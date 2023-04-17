#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS_C 3
#define COLS_O 4
#define C 42

/**
 * Evaluates the performance of TPC-H Q13 (baseline).
 **/

int main(int argc, char** argv) {

  if (argc < 4) {
    printf("\n\nUsage: %s <NUM_ROWS_CUSTOMER> <NUM_ROWS_ORDERS> <BATCH_SIZE>\n\n", argv[0]);
    return -1;
  }

  // initialize communication
  init(argc, argv);

  const int ROWS_C = atol(argv[1]); // input1 size
  const int ROWS_O = atol(argv[2]); // input2 size
  const int BATCH_SIZE_ = atol(argv[3]); // batch size

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS_C, 2*COLS_C, 1};
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {-1, rank, ROWS_O, 2*COLS_O, 2};
  allocate_bool_shares_table(&t2);

  // shares of constant
  BShare s_comm1, s_comm2;

  if (rank == 0) { //P1
    // Initialize input data and shares
    Table r1, r2;
    BShare s_comm3;

    generate_random_table(&r1, ROWS_C, COLS_C);
    generate_random_table(&r2, ROWS_O, COLS_O);
    // generate shares for constant
    generate_bool_share(C, &s_comm1, &s_comm2, &s_comm3);

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS_C, 2*COLS_C, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS_C, 2*COLS_C, 1};
    allocate_bool_shares_table(&t13);

    // t2 Bshare tables for P2, P3 (local to P1)
    BShareTable t22 = {-1, 1, ROWS_O, 2*COLS_O, 2};
    allocate_bool_shares_table(&t22);
    BShareTable t23 = {-1, 2, ROWS_O, 2*COLS_O, 2};
    allocate_bool_shares_table(&t23);

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r1, &t1, &t12, &t13);
    // Generate boolean shares for r2
    generate_bool_share_tables(&r2, &t2, &t22, &t23);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS_C * 2 * COLS_C, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t22.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&s_comm2, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&s_comm3, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS_C * 2 * COLS_C, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t23.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&s_comm3, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&s_comm1, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1.content);
    free(t12.content);
    free(t13.content);
    free(r2.content);
    free(t22.content);
    free(t23.content);

  }
  else { //P2 or P3
    MPI_Recv(&(t1.content[0][0]), ROWS_C * 2 * COLS_C, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&s_comm1, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&s_comm2, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  // start timer
  gettimeofday(&begin, 0);

   // STEP 1: SELECTION ON ORDERS
  #if DEBUG
    if (rank==0) {
      printf("1st SELECTION.\n");
    }
  #endif

  BShare *comm = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(comm!=NULL);
  BShare *rem_comm = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_comm!=NULL);

  // populate vector
  for(int i=0; i<ROWS_O; i++) {
    comm[i] = t2.content[i][4];
    rem_comm[i] = t2.content[i][5];
  }

  BShare *sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(sel!=NULL);
  BShare *rem_sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_sel!=NULL);

  eq_b_array_const(comm, rem_comm, s_comm1, s_comm2, ROWS_O, sel);

  // compute not selected
  for (int i=0; i<ROWS_O; i++) {
      sel[i] ^= 1;
  }
  exchange_shares_array(sel, rem_sel, ROWS_O);

  // Copy arithmetic selection bits to columns 6, 7
  for (int i=0; i<ROWS_O; i++) {
    t2.content[i][6] = sel[i];
    t2.content[i][7] = rem_sel[i];
  }
  free(sel); free(rem_sel); free(comm); free(rem_comm);

  // STEP 2: Join CUSTOMER and ORDERS on CUSTOMERKEY (0, 2)
  BShare *join_res = (BShare *) malloc(ROWS_C*ROWS_O*sizeof(BShare));
  assert(join_res!=NULL);
  BShare *rem_join_res = (BShare *) malloc(ROWS_C*ROWS_O*sizeof(BShare));
  assert(rem_join_res!=NULL);

  Predicate_B p_join = {EQ, NULL, NULL, 0, 2};
  join_b_batch(&t1, &t2, 0, ROWS_C, 0, ROWS_O, p_join, rem_join_res, join_res);
  exchange_shares_array(join_res, rem_join_res, ROWS_C*ROWS_O);

  // STEP 3: AND join result with selection result
  for (int i=0; i<ROWS_C*ROWS_O; i++) {
    join_res[i] = and_b(join_res[i], rem_join_res[i],
                        t2.content[i % ROWS_O][6], t2.content[i % ROWS_O][7],
                        get_next_rb());
  }
  exchange_shares_array(join_res, rem_join_res, ROWS_C*ROWS_O);

  // allocate join result table
  // 0: C_CUSTOMERKEY, 2: O_CUSTOMERKEY, 4: SELECTED (c_count)
  BShareTable res_table = {-1, rank, ROWS_C*ROWS_O, 2*3, 1};
  allocate_bool_shares_table(&res_table);
  for (int i=0; i<ROWS_C; i++) {
    for (int j=0; j<ROWS_O; j++) {
      res_table.content[j + i * ROWS_O][0] = t1.content[i][0];
      res_table.content[j + i * ROWS_O][1] = t1.content[i][1];
      res_table.content[j + i * ROWS_O][2] = t2.content[j][2];
      res_table.content[j + i * ROWS_O][3] = t2.content[j][3];
      res_table.content[j + i * ROWS_O][4] = join_res[j + i * ROWS_O];
      res_table.content[j + i * ROWS_O][5] = rem_join_res[j + i * ROWS_O];
    }
  }
  
  free(join_res); free(rem_join_res);
  // we don't need original tables anymore
  free(t1.content); free(t2.content);

  // STEP 4: SORT on C_CUSTOMERKEY (0)
  #if DEBUG
    if (rank==0) {
      printf("ORDER-BY C_CUSTOMERKEY.\n");
    }
  #endif
  unsigned int att_index[1] = {0};
  bool asc[1] = {1};
  long res_len = ROWS_C*ROWS_O;
  bitonic_sort_batch(&res_table, att_index, 1, asc, res_len/2);

  // STEP 5: GROUP-BY-COUNT on C_CUSTOMERKEY (res_table.0)
  #if DEBUG
    if (rank==0) {
      printf("GROUP-BY on C_CUSTOMERKEY.\n");
    }
  #endif

  unsigned key_indices[1] = {0};
  group_by_sum_rca_odd_even(&res_table, BATCH_SIZE_, key_indices, 1);

  // STEP 6: SORT on C_COUNT (4)
  #if DEBUG
    if (rank==0) {
      printf("ORDER-BY C_COUNT.\n");
    }
  #endif
  att_index[0] = 4;
  bitonic_sort_batch(&res_table, att_index, 1, asc, res_len/2);

  /** Copy C_COUNT to cols 2, 3
   * and initialize last 2 columns to bool shares of 1 for the custdist **/
  for (int i=0; i<res_len; i++) {
    res_table.content[i][2] = res_table.content[i][4];
    res_table.content[i][3] = res_table.content[i][5];
    res_table.content[i][4] = rank % 2;
    res_table.content[i][5] = succ % 2;
  }

  // STEP 7: GROUP-BY-COUNT on C_COUNT (res_table.2)
  #if DEBUG
    if (rank==0) {
      printf("GROUP-BY on C_COUNT.\n");
    }
  #endif

  key_indices[0] = 2;
  group_by_sum_rca_odd_even(&res_table, BATCH_SIZE_, key_indices, 1);

  // STEP 8: FINAL ORDER-BY custdist desc, c_count desc
  #if DEBUG
    if (rank==0) {
      printf("ORDER-BY CUSTDIST, C_COUNT.\n");
    }
  #endif
  unsigned int ind[2] = {4, 2};
  bool asc2[2] = {0, 0};
  bitonic_sort_batch(&res_table, ind, 2, asc2, res_len/2);

  // Open result
  BShare *s_result = (BShare *) malloc(2*res_len*sizeof(BShare));
  Data *result = (Data *) malloc(2*res_len*sizeof(Data));
  for (int i=0; i<res_len; i+=2) {
    s_result[i] = res_table.content[i][2]; // C_COUNT
    s_result[i+1] = res_table.content[i][4]; // CUSTDIST
  }
  open_b_array(s_result, 2*res_len, result);

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("\tTPCH-Q13-BASELINE\t%d\t%d\t%.3f\n", ROWS_C, ROWS_O, elapsed);
  }

  free(res_table.content); free(result); free(s_result);

  // tear down communication
  MPI_Finalize();
  return 0;
}
