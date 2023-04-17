#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193

int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  BShare xs1[10], xs2[10], xs3[10], ys1[10], ys2[10], ys3[10];

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data x[10] = {111, -4, -17, 2345, 999, 0, -28922, 1231241, 0, 23437};
    // Res: 0, 0, 0, 0, 1, 0, 1, 0, 0, 1
    Data y[10] = {0, -4, -5, 123556, 999, 70, -243242, 12421421413421, 0, 78};

    init_sharing();
    for (int i=0; i<10; i++) {
        generate_bool_share(x[i], &xs1[i], &xs2[i], &xs3[i]);
        generate_bool_share(y[i], &ys1[i], &ys2[i], &ys3[i]);
    }
    //Send shares to P2
    MPI_Send(&xs2, 10, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ys2, 10, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&xs3, 10, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ys3, 10, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&xs3, 10, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ys3, 10, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&xs1, 10, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ys1, 10, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&xs1, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ys1, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&xs2, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ys2, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&xs1, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ys1, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&xs2, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ys2, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  // test x[i] > y[i]
  BitShare g, g2;
  for (int i=0; i<10; i++) {
    // test x>y
    g = greater(xs1[i], xs2[i], ys1[i], ys2[i]);
    // reveal the result
    Data out1 = open_b(to_bshare(g));

    g2 = geq(xs1[i], xs2[i], ys1[i], ys2[i]);
    // reveal the result
    Data out2 = open_b(to_bshare(g2));

    // assert and print result
    if (rank == 0) {
      #if DEBUG
        printf("[%d] Greater result (open): %lld\n", rank, out1);
        printf("[%d] Geq result (open): %lld\n", rank, out2);
      #endif
      if (i == 0 || i==6 || i==9) {
        assert(out1 == 1);
      } else {
          assert(out1 == 0);
      }
      if (i == 0 || i==1 || i==4 || i==6 || i==8 || i==9) {
        assert(out2 == 1);
      }
      else {
        assert(out2 == 0);
      }
    }
  }

  // test greater_batch
  BitShare res_batch[10], res_batch2[10]; bool out[10], out3[10];
  greater_batch(xs1, xs2, ys1, ys2, 10, res_batch);
  open_bit_array(res_batch, 10, out);
  geq_batch(xs1, xs2, ys1, ys2, 10, res_batch2);
  open_bit_array(res_batch2, 10, out3);

  if (rank == 0) {
    for (int i=0; i<10; i++) {
      #if DEBUG
        printf("[%d] Batch Greater result (open): %lld\n", rank, out);
        printf("[%d] Batch Geq result (open): %lld\n", rank, out3);
      #endif
      if (i == 0 || i==6 || i==9) {
        assert(out[i] == 1);
      } else {
          assert(out[i] == 0);
      }
      if (i == 0 || i==1 || i==4 || i==6 || i==8 || i==9) {
        assert(out3[i] == 1);
      }
      else {
        assert(out3[i] == 0);
      }
    }
  }

  if (rank==0) {
    printf("TEST GREATER(): OK.\n");
    printf("TEST GEQ(): OK.\n");
    printf("TEST GEQ_BATCH(): OK.\n");
    printf("TEST GREATER_BATCH(): OK.\n");
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
