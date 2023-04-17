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

  BShare x1s1[7], x1s2[7], x1s3[7], x2s1[7], x2s2[7], x2s3[7],
         x3s1[7], x3s2[7], x3s3[7], x4s1[6], x4s2[6], x4s3[6],
         x5s1[6], x5s2[6], x5s3[6], x6s1[6], x6s2[6], x6s3[6],
         x7s1[8], x7s2[8], x7s3[8], x8s1[15], x8s2[15], x8s3[15];

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data x1[7] = {1, 0, 0, 0, 1, 1, 0};
    Data x2[7] = {1, 1, 1, 1, 1, 1, 1};
    Data x3[7] = {0, 0, 0, 0, 0, 0, 0};
    Data x4[6] = {0, 0, 0, 0, 0, 0};
    Data x5[6] = {1, 1, 1, 1, 1, 1};
    Data x6[6] = {1, 0, 0, 0, 1, 1};
    Data x7[8] = {1, 0, 0, 0, 1, 1, 0, 1};
    Data x8[15] = {1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0, 0, 1, 1, 1};

    init_sharing();
    for (int i=0; i<7; i++) {
        generate_bool_share(x1[i], &x1s1[i], &x1s2[i], &x1s3[i]);
        generate_bool_share(x2[i], &x2s1[i], &x2s2[i], &x2s3[i]);
        generate_bool_share(x3[i], &x3s1[i], &x3s2[i], &x3s3[i]);
    }
    for (int i=0; i<6; i++) {
        generate_bool_share(x4[i], &x4s1[i], &x4s2[i], &x4s3[i]);
        generate_bool_share(x5[i], &x5s1[i], &x5s2[i], &x5s3[i]);
        generate_bool_share(x6[i], &x6s1[i], &x6s2[i], &x6s3[i]);
    }
    for (int i=0; i<8; i++) {
        generate_bool_share(x7[i], &x7s1[i], &x7s2[i], &x7s3[i]);
    }
    for (int i=0; i<15; i++) {
        generate_bool_share(x8[i], &x8s1[i], &x8s2[i], &x8s3[i]);
    }
    //Send shares to P2
    MPI_Send(&x1s2, 7, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x2s2, 7, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x3s2, 7, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x4s2, 6, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x5s2, 6, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x6s2, 6, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x7s2, 8, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x8s2, 15, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x1s3, 7, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x2s3, 7, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x3s3, 7, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x4s3, 6, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x5s3, 6, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x6s3, 6, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x7s3, 8, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x8s3, 15, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&x1s3, 7, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x2s3, 7, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x3s3, 7, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x4s3, 6, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x5s3, 6, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x6s3, 6, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x7s3, 8, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x8s3, 15, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x1s1, 7, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x2s1, 7, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x3s1, 7, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x4s1, 6, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x5s1, 6, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x6s1, 6, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x7s1, 8, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&x8s1, 15, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else { //P2 and P3
    MPI_Recv(&x1s1, 7, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x2s1, 7, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x3s1, 7, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x4s1, 6, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x5s1, 6, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x6s1, 6, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x7s1, 8, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x8s1, 15, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x1s2, 7, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x2s2, 7, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x3s2, 7, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x4s2, 6, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x5s2, 6, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x6s2, 6, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x7s2, 8, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&x8s2, 15, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  Data res;
  // Test in-bulk ANDs
  and_b_all(x1s1, x1s2, 7);
  res = open_b(x1s1[0]);
  #if DEBUG
  if (rank==0) {
    printf("[%d] Result 1 (open): %lld\n", rank, res);
  }
  #endif
  if (rank==0) {
    assert(res==0);
  }
  and_b_all(x2s1, x2s2, 7);
  res = open_b(x2s1[0]);
  #if DEBUG
  if (rank==0) {
    printf("[%d] Result 2 (open): %lld\n", rank, res);
  }
  #endif
  if (rank==0) {
    assert(res==1);
  }
  and_b_all(x3s1, x3s2, 7);
  res = open_b(x3s1[0]);
  #if DEBUG
  if (rank==0) {
  printf("[%d] Result 3 (open): %lld\n", rank, res);
  }
  #endif
  if (rank==0) {
    assert(res==0);
  }
  and_b_all(x4s1, x4s2, 6);
  res = open_b(x4s1[0]);
  #if DEBUG
  if (rank==0) {
  printf("[%d] Result 4 (open): %lld\n", rank, res);
  }
  #endif
  if (rank==0) {
    assert(res==0);
  }
  and_b_all(x5s1, x5s2, 6);
  res = open_b(x5s1[0]);
  #if DEBUG
  if (rank==0) {
  printf("[%d] Result 5 (open): %lld\n", rank, res);
  }
  #endif
  if (rank==0) {
    assert(res==1);
  }
  and_b_all(x6s1, x6s2, 6);
  res = open_b(x6s1[0]);
  #if DEBUG
  if (rank==0) {
  printf("[%d] Result 6 (open): %lld\n", rank, res);
  }
  #endif
  if (rank==0) {
    assert(res==0);
  }

  if (rank==0) {
    printf("TEST AND_B_ALL(): OK.\n");
  }
  // Test in bulk ANDs per group
  and_b_all_group(x7s1, x7s2, 4, 2);
  Data out[4];
  open_b_array(x7s1, 8, out);
  #if DEBUG
  if (rank==0) {
  for (int i=0; i<4; i++) {
    printf("[%d] Result 7 (open): %lld\n", rank, out[i]);
  }
  }
  #endif
  if (rank==0) {
    assert(out[0]==0);
    assert(out[1]==0);
    assert(out[2]==1);
    assert(out[3]==0);
  }
  and_b_all_group(x8s1, x8s2, 5, 3);
  Data out2[15];
  open_b_array(x8s1, 15, out2);
  #if DEBUG
  if (rank==0) {
  for (int i=0; i<5; i++) {
    printf("[%d] Result 8 (open): %lld\n", rank, out2[i]);
  }
  }
  #endif
  if (rank==0) {
    assert(out2[0]==0);
    assert(out2[1]==0);
    assert(out2[2]==0);
    assert(out2[3]==0);
    assert(out2[4]==1);
  }

  if (rank==0) {
    printf("TEST AND_B_ALL_GROUP(): OK.\n");
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
