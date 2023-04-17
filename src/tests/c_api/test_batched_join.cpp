#include <stdio.h>
#include <assert.h>
#include <sys/random.h>

#include "../test-utils.h"

#define SHARE_TAG 193
#define ROWS1 1000
#define COLS1 2
#define ROWS2 1000
#define COLS2 2
#define NUM_BATCHES_1 100
#define NUM_BATCHES_2 50

int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();
  
  // r1: first relation, r2: second relation
  BShare *r1s1, *r1s2, *r2s1, *r2s2;
  
  r1s1 = (BShare *) malloc(ROWS1*COLS1*sizeof(BShare));
  r1s2 = (BShare *) malloc(ROWS1*COLS1*sizeof(BShare));
  r2s1 = (BShare *) malloc(ROWS1*COLS1*sizeof(BShare));
  r2s2 = (BShare *) malloc(ROWS1*COLS1*sizeof(BShare));

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data r1[ROWS1][COLS1];
    Data r2[ROWS2][COLS2];
    BShare r1s3[ROWS1][COLS1], r2s3[ROWS2][COLS2];

    // generate random data for r1
    for (int i=0; i<ROWS1; i++) {
      for (int j=0; j<COLS1; j++) {
        r1[i][j] = random();
      }
    }

    // generate random data for r2
    for (int i=0; i<ROWS2; i++) {
      for (int j=0; j<COLS2; j++) {
        r2[i][j] = random();
      }
    }

    printf("Done with initialization.\n");

    init_sharing();

    // generate r1 shares
    for (int i=0; i<ROWS1; i++) {
        for (int j=0; j<COLS1; j++) {
            generate_bool_share(r1[i][j], &r1s1[i*COLS1+j], &r1s2[i*COLS1+j], &r1s3[i][j]);
        }
    }

    // generate r2 shares
    for (int i=0; i<ROWS2; i++) {
        for (int j=0; j<COLS2; j++) {
            generate_bool_share(r2[i][j], &r2s1[i*COLS2+j], &r2s2[i*COLS2+j], &r2s3[i][j]);
        }
    }

    printf("Done with share generation.\n");

    //Send shares to P2
    MPI_Send(r1s2, ROWS1*COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(r2s2, ROWS2*COLS2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r1s3[0][0], ROWS1*COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s3[0][0], ROWS2*COLS2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&r1s3[0][0], ROWS1*COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s3[0][0], ROWS2*COLS2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(r1s1, ROWS1*COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(r2s1, ROWS2*COLS2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else if (rank == 1) { //P2
    MPI_Recv(r1s1, ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r2s1, ROWS2*COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r1s2, ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r2s2, ROWS2*COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(r1s1, ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r2s1, ROWS2*COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r1s2, ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(r2s2, ROWS2*COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  printf("Rank [%d]: seeds exchnaged.\n", rank);

  // allocate BShare tables 
  BShareTable in1 = {-1, rank, ROWS1, 2*COLS1, 1}, in2 = {-1, rank, ROWS2, 2*COLS2, 2};
  allocate_bool_shares_table(&in1);
  allocate_bool_shares_table(&in2);

  // copy shares into the BShareTables
  for (int i=0; i<ROWS1; i++) {
    for (int j=0; j<COLS1; j++) {
        in1.content[i][2 * j] = r1s1[i * COLS1 + j];
        in1.content[i][2 * j + 1] = r1s2[i * COLS1 + j];
    }
  }
  // relation 2
  for (int i=0; i<ROWS2; i++) {
    for (int j=0; j<COLS2; j++) {
        in2.content[i][2 * j] = r2s1[i * COLS2 + j];
        in2.content[i][2 * j + 1] = r2s2[i * COLS2 + j];
    }
  }

  // free temp share tables
  free(r1s1);
  free(r1s2);
  free(r2s1);
  free(r2s2);

  printf("Rank [%d]: tables allocated... Starting computation.\n", rank);
  Predicate_B p = {EQ, NULL, NULL, 0, 0};

  // test batch join
  int batch_size1 = ROWS1 / NUM_BATCHES_1;
  int batch_size2 = ROWS2 / NUM_BATCHES_2;
  BShare res_batched_local[batch_size1*batch_size2]; // batched join result
  BShare res_batched_remote[batch_size1*batch_size2]; // batched join result
  Data out[batch_size1*batch_size2]; // open result

//  int call_ind = 0;

  for (int i=0; i<NUM_BATCHES_1; i++) {
    for (int j=0; j<NUM_BATCHES_2; j++) {
      join_b_batch(&in1, &in2, i*batch_size1, (i+1)*batch_size1,
                    j*batch_size2, (j+1)*batch_size2, p, res_batched_remote, res_batched_local);
      
      open_b_array(res_batched_local, batch_size1*batch_size2, out);
//      if (rank == 0)
//        printf("Done with join batch #%d.\n", ++call_ind);
    }
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
