#include <stdio.h>

#include "mpi.h"
#include "../src/utils.h"
#include "../src/relational.h"

#define SHARE_TAG 193
#define COLS 2

/**
 *
 * SELECT DISTINCT(t1.0)
 * FROM t1, t2
 * WHERE t1.0=t2.0
 *
**/

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("\n\nUsage: %s <NUM_ROWS_1> <NUM_ROWS_2>\n\n", argv[0]);
    return -1;
  }

  // Initialize communication
  init(argc, argv);

  // The party's id
  const int rank = get_rank();

  if (rank == 0) {
    printf("Initializing environment...\n");
  }

  const long ROWS1 = atol(argv[1]); // input1 size
  const long ROWS2 = atol(argv[2]); // input2 size

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS1, 2*COLS, 1};
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {-1, rank, ROWS2, 2*COLS, 2};
  allocate_bool_shares_table(&t2);

  init_tables(&t1, &t2);

  if (rank == 0) {
    printf("Done.\nRunning query...\n");
  }

  // The code below implements the optimization described in Section 4.2
  // of the paper 'Efficient oblivious queries on secret-shared data'
  BShare *join_res = malloc(ROWS1*sizeof(BShare));
  BShare *join_remote = malloc(ROWS1*sizeof(BShare));

  // STEP 1: Sort left input on distinct attribute
  unsigned int sort_att[1] = {0};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, sort_att, 1, asc, t1.numRows/2);

  // STEP 2: Distinct
  BitShare *res_distinct = malloc(ROWS1*sizeof(BitShare));
  distinct_batch(&t1, sort_att[0], res_distinct, t1.numRows - 1);
  // Exchange distinct result shares
  BitShare *res_distinct_remote = malloc(ROWS1*sizeof(BitShare));
  exchange_bit_shares_array(res_distinct, res_distinct_remote, ROWS1);

  // STEP 3: Semi-join
  in(&t1, &t2, 0, 0, join_res, t1.numRows);
  // Exchange semi-join result shares
  exchange_shares_array(join_res, join_remote, ROWS1);

  // STEP 4: Mask tuples before opening
  BShare max=0xFFFFFFFFFFFFFFFF;
  BShare *b_open = malloc(ROWS1*sizeof(BShare));
  assert(b_open !=NULL);
  BShare *b_open_remote = malloc(ROWS1*sizeof(BShare));
  assert(b_open_remote !=NULL);
  BShare *res = malloc(ROWS1*sizeof(BShare));
  assert(res !=NULL);
  for (int i=0; i< ROWS1; i++) { // We only need to scan the outer relation
    b_open[i] = and_b(res_distinct[i], res_distinct_remote[i],
                          join_res[i], join_remote[i],
                          get_next_rb()) & 1;
  }
  exchange_shares_array(b_open, b_open_remote, ROWS1);
  // Multiplex
  for (int i=0; i< ROWS1; i++) {
    BShare b1 = -b_open[i];
    BShare b2 = -b_open_remote[i];
    // res = b_open * att + (1-b) * dummy
    res[i] = and_b(b1, b2,
                      t1.contents[i][0], t1.contents[i][1],
                      get_next_rb());
    res[i] ^= and_b(~b1, ~b2,
                      max, max,
                      get_next_rb());
  }
  free(b_open); free(b_open_remote);

  // STEP 5: Open result
  Data *open_res = malloc(ROWS1*sizeof(Data));
  assert(open_res !=NULL);
  open_b_array(res, ROWS1, open_res);
  // Party 1 prints the result to standard output
  if (rank == 0) {
    printf("Done.\nResult:\n");
    for (int i=0; i<ROWS1; i++) {
      printf("[%d] (bit) = %lld\n", i, open_res[i]);
    }
  }
  // Free memory
  free(join_res); free(join_remote); free(res); free(open_res);
  free(res_distinct); free(res_distinct_remote); free(t1.contents);

  // Tear down communication
  MPI_Finalize();
  return 0;
}
