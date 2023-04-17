#include <stdio.h>
#include <assert.h>

#include "mpi.h"

#include "../test-utils.h"
#include "limits.h"

#define DEBUG 0
#define SHARE_TAG 193

int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  BShare xs1[10], xs2[10], xs3[10], ys1[10], ys2[10], ys3[10];
  
  BShare zs1[20][1], zs2[20][1], zs3[20][1];

  int shift = (sizeof(Data) * 8) - 2;

  Data z[20][1] = {{15}, {20}, {-4}, {1}, {567778}, {-567778},
                  {-120}, {-120}, {0}, {693}, {0}, {0}, {-(1<<shift)}, {-(1<<shift)},
                  {1<<shift}, {1<<shift}, {LONG_LONG_MAX}, {1}, {LONG_LONG_MIN}, {-1}};

  // Initialize input data and shares
  if (rank == 0) { //P1
    // Initialize input data and shares. Last four elements should overflow
    Data x[10] = {15, -4, 567778, -120, 0, 0, -(1<<shift), 1<<shift,
                  LONG_LONG_MAX, LONG_LONG_MIN};
    Data y[10] = {20, 1, -567778, -120, 693, 0, -(1<<shift), 1<<shift, 1, -1};
    
    init_sharing();

    // generate z shares
    for (int i=0; i<20; i++) {
        for (int j=0; j<1; j++) {
            generate_bool_share(z[i][j], &zs1[i][j], &zs2[i][j], &zs3[i][j]);
        }
    }

    for (int i=0; i<10; i++) {
        generate_bool_share(x[i], &xs1[i], &xs2[i], &xs3[i]);
        generate_bool_share(y[i], &ys1[i], &ys2[i], &ys3[i]);
    }
    //Send shares to P2
    MPI_Send(&xs2, 10, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ys2, 10, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&xs3, 10, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ys3, 10, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs2[0][0], 20, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs3[0][0], 20, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&xs3, 10, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ys3, 10, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&xs1, 10, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ys1, 10, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs3[0][0], 20, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs1[0][0], 20, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&xs1, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ys1, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&xs2, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ys2, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&zs1[0][0], 20, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&zs2[0][0], 20, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&xs1, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ys1, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&xs2, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ys2, 10, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&zs1[0][0], 20, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&zs2[0][0], 20, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  // test x[i] == y[i]
  Data res[10] = {35, -3, 0, -240, 693, 0, 0, 0, LONG_LONG_MAX+1,
                    LONG_LONG_MIN-1};
  int len = 187;  // Each boolean addition requires 187 logical ANDs
  BShare rnums[len];
  for (int i=0; i<10; i++) {
    // populate rnums array
    // TODO: We can instead use get_next_rb() in eq_b()
    for (int j=0; j<len; j++) {
      rnums[j] = get_next_rb();
    }
    // test x==y
    BShare first;
    BShare second;
    boolean_addition(xs1[i], xs2[i], ys1[i], ys2[i], &first, &second, rnums);
    // reveal the result
    Data out = open_b(first);

    // assert and print result
    if (rank == 0) {
      #if DEBUG
      printf("[%d] %d Result (open): %lld == %lld\n", rank, i, out, res[i]);
      #endif
      assert(out == res[i]);
    }
  }

  if (rank == 0) {
    printf("TEST RIPPLE CARRY ADDER: OK.\n");
  }

  // test boolean_addition_batch
  BShare res_batched[10];
  Data res_out[10];
  boolean_addition_batch(xs1, xs2, ys1, ys2, res_batched, 10);
  open_b_array(res_batched, 10, res_out);
  if (rank == 0) {
    for (int i=0; i<10; i++) {
      #if DEBUG
      printf("[%d] %d Result-batched (open): %lld == %lld\n", rank, i, res_out[i], res[i]);
      #endif
      assert(res_out[i] == res[i]);
    }
  }

  if (rank == 0) {
    printf("TEST RIPPLE CARRY ADDER (BATCH): OK.\n");
  }

  // test boolean_addition_batch2
  BShare res_batched2[19];
  Data res_out2[19];
  BShareTable t = {-1, rank, 20, 2*1, 1};
  allocate_bool_shares_table(&t);
  for (int i=0; i<20; i++) {
    t.content[i][0] = zs1[i][0];
    t.content[i][1] = zs2[i][0];
  }
  boolean_addition_batch2(t.content, res_batched2, 0, 1, 19, 0);
  open_b_array(res_batched2, 19, res_out2);
  if (rank == 0) {
    for (int i=0; i<19; i++) {
      if (i%2==0) {
        #if DEBUG
        printf("[%d] %d Result-batched (open): %lld == %lld\n", rank, i, res_out2[i], res[i/2]);
        #endif
        assert(res_out2[i] == res[i/2]);
      }
    }
  }

  if (rank == 0) {
    printf("TEST RIPPLE CARRY ADDER (BATCH2): OK.\n");
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
