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
    Data x[10] = {999, -4, 567778, 123556, -754621, 0, 65487, -12885444, 7, 693};
    // First 5 elements are equal to the correspnding x
    Data y[10] = {999, -4, 567778, 123556, -754621, 70, 28922, 45, -45, 0};

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

  // 1. Test element-based equality: x[i] == y[i]
  BShare res;
  for (int i=0; i<10; i++) {

    // test x==y
    res = eq_b(xs1[i], xs2[i], ys1[i], ys2[i]);
    // reveal the result
    Data out = open_b(res);

    // assert and print result
    if (rank == 0) {
      if (i < 5) {
          // equal
          assert(out == 1);
      } else {
          assert(out == 0);
      }
      #if DEBUG
      printf("[%d] Result (open): %lld\n", rank, out);
      #endif
    }
  }

  if (rank==0) {
    printf("TEST EQ_B(): OK.\n");
  }

  // 2. Test array-based equality
  int array_len = 10;
  BShare res_array[array_len];

  eq_b_array(xs1, xs2, ys1, ys2, array_len, res_array);

  Data out_array[array_len];
  open_b_array(res_array, array_len, out_array);

  // assert and print result
  if (rank == 0) {
    for (int i=0; i<array_len; i++) {
      if (i < 5) {
        assert(out_array[i] == 1);
      } else {
        assert(out_array[i] == 0);
      }
      #if DEBUG
      printf("[%d] Array result (open): %lld\n", rank, out_array[i]);
      #endif
    }
    printf("TEST EQ_B_ARRAY(): OK.\n");
  }


  // tear down communication
  MPI_Finalize();
  return 0;
}
