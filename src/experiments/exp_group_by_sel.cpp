#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define SHARE_TAG 193
#define PRIVATE static
#define COLS 1

/**
 * Evaluates the performance of the batched group-by-count operator.
 **/

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("\n\nUsage: %s <NUM_ROWS> <BATCH_SIZE>\n\n", argv[0]);
    return -1;
  }

  // initialize communication
  init(argc, argv);

  const int ROWS = atol(argv[1]); // input size
  const int BATCH_SIZE_ = atol(argv[2]); // batch size

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS, 2*COLS, 1};
  allocate_bool_shares_table(&t1);

  if (rank == 0) { //P1
    // Initialize input data and shares
    Table r1;
    generate_random_table(&r1, ROWS, COLS);

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS, 2*COLS, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS, 2*COLS, 1};
    allocate_bool_shares_table(&t13);

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r1, &t1, &t12, &t13);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS * 2 * COLS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS * 2 * COLS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1.content);
    free(t12.content);
    free(t13.content);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&(t1.content[0][0]), ROWS * 2 * COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&(t1.content[0][0]), ROWS * 2 * COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  AShare *sel_b = (AShare *) malloc(ROWS*sizeof(BShare));
  AShare *remote_sel_b = (AShare *) malloc(ROWS*sizeof(BShare));
  assert(sel_b!=NULL); assert(remote_sel_b!=NULL);

  AShare *counters = (AShare *) malloc(ROWS*sizeof(AShare));
  AShare *remote_counters = (AShare *) malloc(ROWS*sizeof(AShare));
  assert(counters!=NULL); assert(remote_counters!=NULL);
  
  // The number of random numbers we need for the b2a conversion
  int num_rands = ROWS*log2(ROWS) + 1;
  BShare *rand_a = (BShare *) malloc(num_rands*sizeof(BShare));
  BShare *rand_b = (BShare *) malloc(num_rands*sizeof(BShare));
  assert(rand_a!=NULL); assert(rand_b!=NULL);

  // initialize selected bits and counters
  for (int i=0; i<ROWS; i++) {
    sel_b[i] = rank % 2;
    remote_sel_b[i] = succ % 2;
    counters[i] = rank % 2;
    remote_counters[i] = succ % 2;
  }

  // initialize rand bits (all equal to 1)
  for (int i=0; i<num_rands; i++) {
    rand_a[i] = (unsigned int) 1;
    rand_b[i] = (unsigned int) 1;
  }

  /* =======================================================
     1. Measure group-by-wan
  ======================================================== */
  // start timer
  gettimeofday(&begin, 0);

  unsigned key_indices[1] = {0};
  group_by_count_sel_odd_even(&t1, key_indices, 1, BATCH_SIZE_, sel_b, counters, rand_b, rand_a);
  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("%d %d\tGROUP-BY-ODD-EVEN\t%.3f\n", ROWS, BATCH_SIZE_, elapsed);
  }

  /* =======================================================
     2. Measure group-by-lan
  ======================================================== */
  // start timer
  gettimeofday(&begin, 0);

  group_by_count(&t1, key_indices, 1, sel_b, counters, rand_b, rand_a);
  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("%d\tGROUP-BY\t%.3f\n", ROWS, elapsed);
  }

  free(t1.content); free(counters); free(remote_counters);
  free(sel_b); free(remote_sel_b);
  free(rand_a); free(rand_b);

  // tear down communication
  MPI_Finalize();
  return 0;
}
