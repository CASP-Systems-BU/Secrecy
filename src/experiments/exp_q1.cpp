#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS1 3
#define COLS2 1

/**
 * Evaluates the performance of Q1 (comorbidity).
 **/

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("\n\nUsage: %s <NUM_ROWS_1> <NUM_ROWS_2>\n\n", argv[0]);
    return -1;
  }

  // initialize communication
  init(argc, argv);

  const int ROWS1 = atol(argv[1]); // input1 size
  const int ROWS2 = atol(argv[2]); // input2 size

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS1, 2*COLS1, 1}; // {pid, diag, cnt}
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {-1, rank, ROWS2, 2*COLS2, 2};
  allocate_bool_shares_table(&t2);

  BShare *rb = (BShare *) malloc(ROWS1*sizeof(BShare));
  AShare *ra = (AShare *) malloc(ROWS1*sizeof(AShare));
  int num_rands = ROWS1*log2(ROWS1) + 1;
  BShare *rand_b = (BShare *) malloc(num_rands*sizeof(BShare));
  assert(rand_b!=NULL);
  AShare *rand_a = (AShare *) malloc(num_rands*sizeof(AShare));
  assert(rand_a!=NULL);

  // initialize rand bits for conversion (all equal to 0)
  for (int i=0; i<ROWS1; i++) {
    ra[i] = (unsigned int) 0;
    rb[i] = (unsigned int) 0;
  }

  // initialize rand bits for group-by (all equal to 0)
  for (int i=0; i<num_rands; i++) {
    rand_a[i] = (unsigned int) 0;
    rand_b[i] = (unsigned int) 0;
  }

  if (rank == 0) { //P1
    // Initialize input data and shares
    Table r1, r2;
    generate_random_table(&r1, ROWS1, COLS1);
    generate_random_table(&r2, ROWS2, COLS2);

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS1, 2*COLS1, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS1, 2*COLS1, 1};
    allocate_bool_shares_table(&t13);

    // t2 Bshare tables for P2, P3 (local to P1)
    BShareTable t22 = {-1, 1, ROWS2, 2*COLS2, 2};
    allocate_bool_shares_table(&t22);
    BShareTable t23 = {-1, 2, ROWS2, 2*COLS2, 2};
    allocate_bool_shares_table(&t23);

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r1, &t1, &t12, &t13);
    // Generate boolean shares for r2
    generate_bool_share_tables(&r2, &t2, &t22, &t23);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t22.content[0][0]), ROWS2 * 2 * COLS2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t23.content[0][0]), ROWS2 * 2 * COLS2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1.content);
    free(t12.content);
    free(t13.content);
    free(r2.content);
    free(t22.content);
    free(t23.content);

  }
  else if (rank == 1) { //P2
    MPI_Recv(&(t1.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS2 * 2 * COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&(t1.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS2 * 2 * COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  // start timer
  gettimeofday(&begin, 0);

  // STEP 1: SORT t1 on diag (att=2)
  #if DEBUG
    if (rank==0) {
      printf("Sorting.\n");
    }
  #endif
  unsigned int att_index[1] = {2};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS1/2);

  // STEP 2: IN
  #if DEBUG
    if (rank==0) {
      printf("IN.\n");
    }
  #endif
  BShare *in_res = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(in_res!=NULL);
  in(&t1, &t2, 0, 0, in_res, ROWS1);

  // STEP 3: GROUP-BY-COUNT on diag (att=2)
  #if DEBUG
    if (rank==0) {
      printf("Conversion.\n");
    }
  #endif
  // a. get arithmetic shares of selected bits
  AShare *in_res_a = (AShare *) malloc(ROWS1*sizeof(AShare));
  assert(in_res_a!=NULL);
  convert_single_bit_array(in_res, ra, rb, ROWS1, in_res_a);
  #if DEBUG
    if (rank==0) {
      printf("Group-by.\n");
    }
  #endif
  unsigned key_indices[1] = {2};
  group_by_count_sel_odd_even(&t1, key_indices, 1, ROWS1, in_res, in_res_a, rand_a, rand_b);
  free(rand_a); free(rand_b);

  // reuse ra, rb arrays for exchange of arithmetic counts
  exchange_a_shares_array(in_res_a, ra, ROWS1);

  // STEP 4: sort group's output on count
  // reuse rb, in_res for result of conversion to binary
  #if DEBUG
    if (rank==0) {
      printf("Conversion.\n");
    }
  #endif
  convert_a_to_b_array(in_res_a, ra, rb, in_res, ROWS1);
  free(in_res_a); free(ra);

  // order by count
  // first copy the binary counter to the last column of t1
  for (int i=0; i<ROWS1; i++) {
    t1.content[i][4] = rb[i]; // local share of count
    t1.content[i][5] = in_res[i]; // remote share of count
  }
  free(rb); free(in_res);

  // STEP 5: Sort by cnt
  #if DEBUG
    if (rank==0) {
      printf("Sorting.\n");
    }
  #endif
  att_index[0] = 4;
  asc[0] = 0;
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS1/2);

  Data result[10][2];

  // Open first 10 elements
  for (int i=0; i<10; i++) {
    result[i][0] = open_b(t1.content[i][2]); // diag
    result[i][1] = open_b(t1.content[i][4]); // count
  }

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("\tQ1\t%d\t%d\t%.3f\n", ROWS1, ROWS2, elapsed);
  }

  #if DEBUG
    if (rank == 0) {
      for (int i=0; i<10; i++) {
        printf("[%d] (diag, cnt) = %lld, %lld\n", i, result[i][0], result[i][1]);
      }
    }
  #endif

  free(t1.content); free(t2.content);

  // tear down communication
  MPI_Finalize();
  return 0;
}
