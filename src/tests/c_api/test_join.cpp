#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define ROWS1 5
#define COLS1 2
#define ROWS2 4
#define COLS2 2

int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // r1: first relation with 5 rows, r2: second relation with 4 rows
  BShare r1s1[ROWS1][COLS1], r1s2[ROWS1][COLS1], r1s3[ROWS1][COLS1],
         r2s1[ROWS2][COLS2], r2s2[ROWS2][COLS2], r2s3[ROWS2][COLS2];

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data r1[ROWS1][COLS1] = {{1, 42}, {2, 42}, {3, 42}, {4, 42}, {5, 42}};
    Data r2[ROWS2][COLS2] = {{1, 99}, {3, 99}, {5, 99}, {7, 99}};

    init_sharing();

    // generate r1 shares
    for (int i=0; i<ROWS1; i++) {
        for (int j=0; j<COLS1; j++) {
            generate_bool_share(r1[i][j], &r1s1[i][j], &r1s2[i][j], &r1s3[i][j]);
        }
    }

    // generate r2 shares
    for (int i=0; i<ROWS2; i++) {
        for (int j=0; j<COLS2; j++) {
            generate_bool_share(r2[i][j], &r2s1[i][j], &r2s2[i][j], &r2s3[i][j]);
        }
    }

    //Send shares to P2
    MPI_Send(&r1s2[0][0], 5*2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s2[0][0], 4*2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r1s3[0][0], 5*2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s3[0][0], 4*2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&r1s3[0][0], 5*2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s3[0][0], 4*2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r1s1[0][0], 5*2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s1[0][0], 4*2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&r1s1[0][0], 5*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r2s1[0][0], 4*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r1s2[0][0], 5*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r2s2[0][0], 4*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&r1s1[0][0], 5*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r2s1[0][0], 4*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r1s2[0][0], 5*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r2s2[0][0], 4*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  BShare res[ROWS1*ROWS2]; // join result

  // allocate BShare tables
  BShareTable in1 = {-1, rank, ROWS1, 2*COLS1, 1}, in2 = {-1, rank, ROWS2, 2*COLS2, 2};
  allocate_bool_shares_table(&in1);
  allocate_bool_shares_table(&in2);

  // copy shares into the BShareTables
  for (int i=0; i<ROWS1; i++) {
      in1.content[i][0] = r1s1[i][0];
      in1.content[i][1] = r1s2[i][0];
      in1.content[i][2] = r1s1[i][1];
      in1.content[i][3] = r1s2[i][1];
  }
  // relation 2
  for (int i=0; i<ROWS2; i++) {
      in2.content[i][0] = r2s1[i][0];
      in2.content[i][1] = r2s2[i][0];
      in2.content[i][2] = r2s1[i][1];
      in2.content[i][3] = r2s2[i][1];
  }

  Predicate_B p = {EQ, NULL, NULL, 0, 0};

  join_b(in1, in2, p, res);

  // test batch join
  BShare res_batched[ROWS1*ROWS2], remote[ROWS1*ROWS2]; // batched join result

  join_b_batch(&in1, &in2, 0, 3, 0, 2, p, remote, res_batched);
  join_b_batch(&in1, &in2, 3, 5, 0, 2, p, &remote[6], &res_batched[6]);

  join_b_batch(&in1, &in2, 0, 3, 2, 4, p, &remote[10], &res_batched[10]);
  join_b_batch(&in1, &in2, 3, 5, 2, 4, p, &remote[16], &res_batched[16]);

  // reveal the result
  Data out[ROWS1*ROWS2]; /** TODO: only rank0 needs to allocate **/
  open_b_array(res, ROWS1*ROWS2, out);

  // assert and print result
  if (rank == 0) {
      for (int i=0; i<ROWS1*ROWS2; i++) {
          #if DEBUG
            printf("[%d] %lld\t", i, out[i]);
          #endif
          if (i==0 || i==9 || i==18) {
              assert(out[i] == 1);
          }
          else {
              assert(out[i] == 0);
          }
      }
      printf("TEST JOIN: OK.\n");
  }

  // reveal the batched result
  open_b_array(res_batched, ROWS1*ROWS2, out);

  // assert and print result
  if (rank == 0) {
      for (int i=0; i<ROWS1*ROWS2; i++) {
          #if DEBUG
            printf("[%d] %lld\t", i, out[i]);
          #endif
          if (i==0 || i==5 || i==18) {
              assert(out[i] == 1);
          }
          else {
              assert(out[i] == 0);
          }
      }
      printf("TEST BATCH JOIN: OK.\n");
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
