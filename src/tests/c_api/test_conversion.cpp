#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define ROWS 10

/**
 * Tests single-bit conversion from binary to arithmetic sharing
**/
int main(int argc, char** argv) {

  // initialize communication and sharing
  init(argc, argv);
  init_sharing();

  const int rank = get_rank();

  Data r[ROWS] = {1, 0, 0, 0, 1, 1, 1, 1, 0, 1};

  AShare converted[ROWS];
  BShare rs[ROWS], rb[ROWS];
  AShare ra[ROWS];

  // P1 generates Data BShares and random bit shares
  if (rank == 0) {

    BShare rs2[ROWS], rs3[ROWS], rb2[ROWS], rb3[ROWS];
    AShare ra2[ROWS], ra3[ROWS];

    for (int i=0; i<ROWS; i++) {
      generate_bool_share(r[i], &rs[i], &rs2[i], &rs3[i]);
    }

    //Send shares to P2
    MPI_Send(&rs2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&rs3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // Generate random bits and corresponding shares
    generate_rand_bit_shares(rb, ra, rb2, ra2, rb3, ra3, ROWS);

    // Send random bit shares
    //Send shares to P2
    MPI_Send(&rb2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ra2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&rb3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ra3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&rs, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&rb, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ra, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&rs, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&rb, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ra, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  for (int i=0; i<ROWS; i++) {
    converted[i] = convert_single_bit(rs[i], ra[i], rb[i]);
  }

  // reveal the result
  Data out[ROWS];
  open_array(converted, ROWS, out);

  // assert and print result
  if (rank == 0) {
    for (int i=0; i<ROWS; i++) {
        #if DEBUG
          printf("[%d] %lld\t", i, out[i]);
        #endif
        if (i==0 || (i > 3 && i < 8) || i==9) {
          assert(out[i] == 1);
        }
        else {
          assert(out[i] == 0);
        }
    }
    printf("TEST B2A CONVERSION (SINGLE-BIT): OK.\n");
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
