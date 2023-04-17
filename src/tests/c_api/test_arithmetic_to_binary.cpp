#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define ROWS 10

/**
 * Tests conversion from arithmetic sharing to binary
**/
int main(int argc, char** argv) {

  // initialize communication and sharing
  init(argc, argv);
  init_sharing();

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  Data r[ROWS] = {10, -34, 2, 34, 778, 2334, 0, -654, 12222, 9};

  BShare converted_local[ROWS], converted_remote[ROWS];
  AShare local[ROWS], remote[ROWS];

  // P1 generates Data BShares and random bit shares
  if (rank == 0) {

    AShare ra2[ROWS], ra3[ROWS];

    for (int i=0; i<ROWS; i++) {
      generate_int_share(r[i], &local[i], &ra2[i], &ra3[i]);
    }

    //Send shares to P2
    MPI_Send(&ra2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&ra3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&local, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&local, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  exchange_a_shares_array(local, remote, ROWS);
  AShare prev_local[ROWS];
  for (int i=0; i<ROWS; i++) {
    prev_local[i] = local[i];
  }
  convert_a_to_b_array(local, remote, converted_local, converted_remote, ROWS);

  // reveal the result
  Data out_a[ROWS];
  Data out_b[ROWS];
  open_array(prev_local, ROWS, out_a);
  open_b_array(converted_local, ROWS, out_b);

  // assert and print result
  if (rank == 0) {
    for (int i=0; i<ROWS; i++) {
      #if DEBUG
        printf("[%d] %lld\n", i, out_b[i]);
      #endif
      assert(out_a[i] == out_b[i]);
    }
    printf("TEST A2B CONVERSION: OK.\n");
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
