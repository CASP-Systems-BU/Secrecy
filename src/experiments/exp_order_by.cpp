#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define SHARE_TAG 193
#define PRIVATE static

/**
 * Evaluates the performance of ORDER-BY (bitonic sort).
 **/

int main(int argc, char** argv) {

  if (argc < 4) {
    printf("\n\nUsage: %s <NUM_ROWS> <NUM_COLS> <NUM_ATTRIBUTES>\n\n", argv[0]);
    return -1;
  }

  // initialize communication
  init(argc, argv);

  const int ROWS = atol(argv[1]); // input size
  const int COLS = atoi(argv[2]); // number of columns
  const int ATT = atoi(argv[3]); // number of sorting attributes (up to 3)
  assert(ATT<=COLS);
  assert(ATT<4 && ATT>0);

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

  unsigned int att_index[3] = {0, 2, 4};
  bool asc[3] = {1, 1, 1};
   
  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  /* =======================================================
     1. Measure CMP-SWAP batched
  ======================================================== */
  // start timer
  gettimeofday(&begin, 0);

  bitonic_sort_batch(&t1, att_index, ATT, asc, ROWS/2);
  
  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("%d\tORDER-BY\t%d\t%.3f\n", ATT, ROWS, elapsed);
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}