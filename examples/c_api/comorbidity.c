#include <stdio.h>

#include "mpi.h"
#include "../src/utils.h"
#include "../src/relational.h"

#define SHARE_TAG 193
#define COLS1 3
#define COLS2 1

/**
 * Comorbidity query.
 *
 * SELECT diag, COUNT(*) as cnt
 * FROM diagnosis
 * WHERE pid IN cdiff_cohort
 * GROUP BY diag
 * ORDER BY cnt DESC
 * LIMIT 10
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
  BShareTable t1 = {-1, rank, ROWS1, 2*COLS1, 1}; // Diagnosis(pid, diag, cnt)
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {-1, rank, ROWS2, 2*COLS2, 2}; // cdiff_cohort(pid)
  allocate_bool_shares_table(&t2);
  // Arrays used for dual sharing
  BShare *rb = calloc(ROWS1, sizeof(BShare));
  AShare *ra = calloc(ROWS1, sizeof(AShare));
  BShare *rand_b = calloc(2*(ROWS1-1), sizeof(BShare));
  AShare *rand_a = calloc(2*(ROWS1-1), sizeof(AShare));

  // Initialize and distribute secret shares
  init_tables(&t1, &t2);

  if (rank == 0) {
    printf("Done.\nRunning query...\n");
  }

  // STEP 1: SORT t1 on diag (att=2)
  unsigned int att_index[1] = {2};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS1/2);

  // STEP 2: IN
  BShare *in_res = malloc(ROWS1*sizeof(BShare));
  assert(in_res!=NULL);
  in(&t1, &t2, 0, 0, in_res, ROWS1);

  // STEP 3: GROUP-BY-COUNT on diag (att=2)
  AShare *in_res_a = malloc(ROWS1*sizeof(AShare));
  assert(in_res_a!=NULL);
  // Dual sharing: Get arithmetic shares of selected bits
  convert_single_bit_array(in_res, ra, rb, ROWS1, in_res_a);
  group_by_count(&t1, 2, in_res, in_res_a, rand_a, rand_b);
  free(rand_a); free(rand_b);
  exchange_a_shares_array(in_res_a, ra, ROWS1);
  convert_a_to_b_array(in_res_a, ra, rb, in_res, ROWS1);
  free(in_res_a); free(ra);

  // Copy the boolean counter to the last column of t1
  for (int i=0; i<ROWS1; i++) {
    t1.contents[i][4] = rb[i]; // local share of count
    t1.contents[i][5] = in_res[i]; // remote share of count
  }
  free(rb); free(in_res);

  // STEP 4: Sort by cnt in descending order
  att_index[0] = 4;
  asc[0] = 0;
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS1/2);

  // STEP 5: Open first 10 rows in the result
  Data result[10][2];
  int limit = 10 <= t1.numRows ? 10 : t1.numRows;
  for (int i=0; i<limit; i++) {
    result[i][0] = open_b(t1.contents[i][2]); // diag
    result[i][1] = open_b(t1.contents[i][4]); // count
  }
  // Party 1 prints the result to standard output
  if (rank == 0) {
    printf("Done.\nResult:\n");
    for (int i=0; i<10; i++) {
      printf("[%d] (diag, cnt) = %lld, %lld\n", i, result[i][0], result[i][1]);
    }
  }
  free(t1.contents); free(t2.contents);

  // Tear down communication
  MPI_Finalize();
  return 0;
}
