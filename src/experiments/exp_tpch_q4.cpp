#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS_O 5
#define COLS_L 4
#define D1 1993
#define D2 1995

/**
 * Evaluates the performance of TPC-H Q4.
 **/

int main(int argc, char** argv) {

  if (argc < 4) {
    printf("\n\nUsage: %s <NUM_ROWS_ORDERS> <NUM_ROWS_LINEITEM> <BATCH_SIZE>\n\n", argv[0]);
    return -1;
  }

  // initialize communication
  init(argc, argv);

  const int ROWS_O = atol(argv[1]); // input1 size
  const int ROWS_L = atol(argv[2]); // input2 size
  const int BATCH_SIZE_ = atol(argv[3]); //batch size for semi-join

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS_O, 2*COLS_O, 1};
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {-1, rank, ROWS_L, 2*COLS_L, 2};
  allocate_bool_shares_table(&t2);

  BShare *rb = (BShare *) calloc(ROWS_O, sizeof(BShare));
  AShare *ra = (AShare *) calloc(ROWS_O, sizeof(AShare));
  int num_rands = ROWS_O*log2(ROWS_O) + 1;
  BShare *rand_b = (BShare *) calloc(num_rands, sizeof(BShare));
  AShare *rand_a = (AShare *) calloc(num_rands, sizeof(AShare));

  if (rank == 0) { //P1
    // Initialize input data and shares
    Table r1, r2;
    generate_random_table(&r1, ROWS_O, COLS_O);
    generate_random_table(&r2, ROWS_L, COLS_L);

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS_O, 2*COLS_O, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS_O, 2*COLS_O, 1};
    allocate_bool_shares_table(&t13);

    // t2 Bshare tables for P2, P3 (local to P1)
    BShareTable t22 = {-1, 1, ROWS_L, 2*COLS_L, 2};
    allocate_bool_shares_table(&t22);
    BShareTable t23 = {-1, 2, ROWS_L, 2*COLS_L, 2};
    allocate_bool_shares_table(&t23);

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r1, &t1, &t12, &t13);
    // Generate boolean shares for r2
    generate_bool_share_tables(&r2, &t2, &t22, &t23);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t22.content[0][0]), ROWS_L * 2 * COLS_L, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t23.content[0][0]), ROWS_L * 2 * COLS_L, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1.content);
    free(t12.content);
    free(t13.content);
    free(r2.content);
    free(t22.content);
    free(t23.content);

  }
  else { //P2 or P3
    MPI_Recv(&(t1.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS_L * 2 * COLS_L, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  // start timer
  gettimeofday(&begin, 0);

  // STEP 1: SORT ORDERS on priority (att=4)
  #if DEBUG
    if (rank==0) {
      printf("Sorting.\n");
    }
  #endif
  unsigned int att_index[1] = {4};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS_O/2);

  // STEP 2: Selection L_COMMITDATE < L_RECEIPTDATE
  #if DEBUG
    if (rank==0) {
      printf("1st selection on LINEITEM.\n");
    }
  #endif

  BShare *sel_l = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(sel_l!=NULL);
  BShare *rem_sel_l = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(rem_sel_l!=NULL);

  Predicate_B p_l = {GT, NULL, NULL, 6, -1};
  select_b(t2, p_l, sel_l);
  exchange_shares_array(sel_l, rem_sel_l, ROWS_L);

  // Copy selection bits to columns 6, 7
  for (int i=0; i<ROWS_L; i++) {
    t2.content[i][6] = sel_l[i];
    t2.content[i][7] = rem_sel_l[i];
  }

  free(sel_l); free(rem_sel_l);

  // STEP 3: Fused Selection-IN
  #if DEBUG
    if (rank==0) {
      printf("(fused) IN.\n");
    }
  #endif

  BShare *in_res = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(in_res!=NULL);
  BShare *rem_in_res = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_in_res!=NULL);

  in_sel_right(&t1, &t2, 0, 0, 6, in_res, BATCH_SIZE_);
  exchange_shares_array(in_res, rem_in_res, ROWS_O);

  // STEP 4: Apply selection O_ORDERDATE >= D1
  #if DEBUG
    if (rank==0) {
      printf("1st selection on ORDERS.\n");
    }
  #endif

  BShare *sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(sel!=NULL);
  BShare *rem_sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_sel!=NULL);

  Predicate_B p = {GEQ, NULL, NULL, 6, -1};
  select_b(t1, p, sel);
  exchange_shares_array(sel, rem_sel, ROWS_O);

  // Copy selection bits to O_ORDERDATE - D1 column
  for (int i=0; i<ROWS_O; i++) {
    t1.content[i][6] = sel[i];
    t1.content[i][7] = rem_sel[i];
  }

  // STEP 5: Apply selection O_ORDERDATE < D2
  #if DEBUG
    if (rank==0) {
      printf("2nd selection on ORDERS.\n");
    }
  #endif

  Predicate_B p2 = {GT, NULL, NULL, 8, -1};
  select_b(t1, p2, sel);
  exchange_shares_array(sel, rem_sel, ROWS_O);

  // Copy selection bits to O_ORDERDATE - D2 column
  for (int i=0; i<ROWS_O; i++) {
    t1.content[i][8] = sel[i];
    t1.content[i][9] = rem_sel[i];
  }

  // STEP 6: Compute AND of selections
  BShare mask = 1;
  #if DEBUG
    if (rank==0) {
      printf("AND of selections.\n");
    }
  #endif

  BShare *res_sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(res_sel!=NULL);
  BShare *rem_res_sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_res_sel!=NULL);

  for (int j=0; j<ROWS_O; j++) {
    res_sel[j] = and_b(t1.content[j][6], t1.content[j][7],
                       t1.content[j][8], t1.content[j][9], get_next_rb())
                      & mask;
  }
  exchange_shares_array(res_sel, rem_res_sel, ROWS_O);

  #if DEBUG
    if (rank==0) {
      printf("AND with result of IN.\n");
    }
  #endif

  for (int j=0; j<ROWS_O; j++) {
    res_sel[j] = and_b(in_res[j], rem_in_res[j],
                      res_sel[j], rem_res_sel[j], get_next_rb())
                      & mask;
  }

  free(in_res); free(rem_in_res);

  // STEP 7: GROUP-BY-COUNT on O_ORDERPRIORITY (att=4)
  #if DEBUG
    if (rank==0) {
      printf("Conversion.\n");
    }
  #endif
  // a. get arithmetic shares of selected bits
  AShare *sel_a = (AShare *) malloc(ROWS_O*sizeof(AShare));
  assert(sel_a!=NULL);
  convert_single_bit_array(res_sel, ra, rb, ROWS_O, sel_a);
  #if DEBUG
    if (rank==0) {
      printf("Group-by.\n");
    }
  #endif
  unsigned key_indices[1] = {4};
  group_by_count_sel_odd_even(&t1, key_indices, 1, BATCH_SIZE_, res_sel, sel_a, rand_b, rand_a);

  free(ra); free(rb); free(sel); free(rem_sel); free(res_sel);
  free(rand_a); free(rand_b);

  // Open counts
  Data *res_count = (Data *) malloc(ROWS_O*sizeof(Data));
  open_array(sel_a, ROWS_O, res_count);
  // Open O_ORDERPRIORITY
  Data *res_prio = (Data *) malloc(ROWS_O*sizeof(Data));
  for (int i=0; i<ROWS_O; i++) {
    sel_a[i] = t1.content[i][4]; // O_ORDERPRIORITY
  }
  open_b_array(sel_a, ROWS_O, res_prio);

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("\tTPCH-Q4\t%d\t%d\t%.3f\n", ROWS_O, ROWS_L, elapsed);
  }

  free(t1.content); free(t2.content); free(sel_a); free(res_count); free(res_prio);

  // tear down communication
  MPI_Finalize();
  return 0;
}
