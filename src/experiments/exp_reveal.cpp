#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

/**
 * Evaluates the performance of reveal_array.
 **/

int main(int argc, char** argv) {

  if (argc < 2) {
    printf("\n\nUsage: %s [INPUT_SIZE]\n\n", argv[0]);
    return -1;
  }

  // initialize communication
  init(argc, argv);

  const int ROWS = atol(argv[1]); // input size

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();
  
  BShare *r1s1, *r1s2;
  r1s1 = (BShare *) malloc(ROWS*sizeof(BShare));
  r1s2 = (BShare *) malloc(ROWS*sizeof(BShare));

  generate_and_share_random_data(rank, r1s1, r1s2, ROWS);

  //exchange seeds
  exchange_rsz_seeds(succ, pred);
   
  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  /* =======================================================
     1. Measure reveal_array synchronous
  ======================================================== */
  // start timer
  gettimeofday(&begin, 0);

  reveal_b_array(r1s1, ROWS);
  
  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("SYNC\t%d\t%.3f\n", ROWS, elapsed);
  }

  /* =======================================================
     2. Measure reveal_array asynchronous
  ======================================================== */
  // start timer
  gettimeofday(&begin, 0);

  reveal_b_array_async(r1s1, ROWS);
  
  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("ASYNC\t%d\t%.3f\n", ROWS, elapsed);
  }

  free(r1s1); free(r1s2);

  // tear down communication
  MPI_Finalize();
  return 0;
}
