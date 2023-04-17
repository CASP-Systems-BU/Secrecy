/**
 *
 * This program tests random sharing of zero among 3 parties.
 *
 * Each party generates a random seed and sends it to its successor:
 * P0 to P1, P1 to P2, P2 to P0.
 *
 * Then, the party can use its local seed and the remote seed
 * to generate two random numbers: r_local, r_remote.
 *
 * The random number to be used for multiplication is:
 * R = r_local - r_remote, so that R0 ^ R1 ^ R2 = 0.
 *
 **/
#include "mpi.h"
#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

void test_next_rb(int size);
void test_next_rb_array(int size);

int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  // exchange random seeds
  exchange_rsz_seeds(get_succ(), get_pred());

  test_next_rb(100);
  test_next_rb_array(500);

  // tear down communication
  MPI_Finalize();
  return 0;
}
// test individual binary r generation
void test_next_rb(int size) {
  BShare rs[size];
  for (int i=0; i<size; i++) {
    rs[i] = get_next_rb();
  }

  // open and have rank 0 test the results
  Data *out = rs;
  open_b_array(rs, size, out);

  if (get_rank() == 0) {
    for (int i=0; i<size; i++) {
      assert(out[i] == 0);
    }
    printf("TEST GET_NEXT_RB(): OK.\n");
  }
}

// test generation of an array of binary rnums
void test_next_rb_array(int size) {
  BShare rs[size];
  get_next_rb_array(rs, size);

  // open and have rank 0 test the results
  Data *out = rs;
  open_b_array(rs, size, out);

  if (get_rank() == 0) {
    for (int i=0; i<size; i++) {
      assert(out[i] == 0);
    }
    printf("TEST GET_NEXT_RB_ARRAY(): OK.\n");
  }
}
