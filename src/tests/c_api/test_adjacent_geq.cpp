#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define ROWS 10

int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  BShare zs1[10][2], zs2[10][2], zs3[10][2];

  // Initialize input data and shares

  Data z[10][2] = {{2, 42}, {1, 42}, {2, 42}, {3, 42}, {15, 42}, {15, 43},
                  {15, 44}, {17, 1}, {19, 1}, {18, 1}};

  if (rank == 0) { //P1

    init_sharing();

    // generate z shares
    for (int i=0; i<10; i++) {
        for (int j=0; j<2; j++) {
            generate_bool_share(z[i][j], &zs1[i][j], &zs2[i][j], &zs3[i][j]);
        }
    }

    //Send shares to P2
    MPI_Send(&zs2[0][0], 10*2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs3[0][0], 10*2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&zs3[0][0], 10*2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs1[0][0], 10*2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&zs1[0][0], 10*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&zs2[0][0], 10*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&zs1[0][0], 10*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&zs2[0][0], 10*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  // test group_by_count
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing adjacent_geq\n", rank);
    }
  #endif

  BShareTable t2 = {-1, rank, 10, 2*2, 1};
  allocate_bool_shares_table(&t2);
  // copy shares into the BShareTables
  for (int i=0; i<10; i++) {
      t2.content[i][0] = zs1[i][0];
      t2.content[i][1] = zs2[i][0];
      t2.content[i][2] = zs1[i][1];
      t2.content[i][3] = zs2[i][1];
  }

  BitShare res[ROWS-1];

  // Check if r_i[0] <= r_{i+1}[0]
  adjacent_geq(&t2, 0, 0, res, 7, 1);

  // reveal the result
  bool out[ROWS-1];
  open_bit_array(res, ROWS-1, out);

  if (rank==0) {
    // Data z[10][2] = {{2, 42}, {1, 42}, {2, 42}, {3, 42}, {15, 42}, {15, 43},
    //                 {15, 44}, {17, 1}, {19, 1}, {18, 1}};
    Data s[ROWS-1] = {0,1,1,1,1,1,1,1,0};

    #if DEBUG
      printf("[%d] Adjacent geq array (open):\n", rank);
    #endif
    for (int i=0; i<ROWS-1; i++) {
      #if DEBUG
        printf("%d \t", out[i]);
      #endif
      assert(out[i]==s[i]);
    }
    printf("TEST ADJACENT GEQ: OK.\n");
  }
  // tear down communication
  MPI_Finalize();
  return 0;
}
