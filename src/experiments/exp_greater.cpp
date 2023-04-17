#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define SHARE_TAG 193

/**
 * Evaluates the performance of greater_batch.
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
  
  BShare *r1s1, *r1s2, *r2s1, *r2s2;
  r1s1 = (BShare *) malloc(ROWS*sizeof(BShare));
  r1s2 = (BShare *) malloc(ROWS*sizeof(BShare));
  r2s1 = (BShare *) malloc(ROWS*sizeof(BShare));
  r2s2 = (BShare *) malloc(ROWS*sizeof(BShare));

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data *r1, *r2;
    r1 = (Data *) malloc(ROWS*sizeof(Data));
    r2 = (Data *) malloc(ROWS*sizeof(Data));
    BShare *r1s3, *r2s3;
    r1s3 = (BShare *) malloc(ROWS*sizeof(BShare));
    r2s3 = (BShare *) malloc(ROWS*sizeof(BShare));

    // generate random data for r1 and r2
    for (long i=0; i<ROWS; i++) {
      r1[i] = random();
      r2[i] = random();
    }

    init_sharing();

    // generate r1 and r2 shares
    for (long i=0; i<ROWS; i++) {
      generate_bool_share(r1[i], &r1s1[i], &r1s2[i], &r1s3[i]);
      generate_bool_share(r2[i], &r2s1[i], &r2s2[i], &r2s3[i]);
    }

    //Send shares to P2
    MPI_Send(r1s2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(r1s3, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(r2s2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(r2s3, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(r1s3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(r1s1, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(r2s3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(r2s1, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1);
    free(r2);
    free(r1s3);
    free(r2s3);
  }
  else if (rank == 1) { //P2
    MPI_Recv(r1s1, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r1s2, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r2s1, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r2s2, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(r1s1, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r1s2, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r2s1, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r2s2, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);
   
  struct timeval begin, end;
  long seconds, micro;
  double elapsed;
  
  /* =======================================================
     1. Measure array-based greater (async) 
  ======================================================== */
  BitShare *res_array = (BitShare *) malloc(ROWS*sizeof(BitShare));

  // start timer
  gettimeofday(&begin, 0);

  greater_batch(r1s1, r1s2, r2s1, r2s2, ROWS, res_array);

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("%d\t%.3f\n", ROWS, elapsed);
  }

  free(r1s1); free(r1s2); free(r2s1); free(r2s2); free(res_array);

  // tear down communication
  MPI_Finalize();
  return 0;
}
