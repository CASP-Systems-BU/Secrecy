#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define SHARE_TAG 193
#define PRIVATE static
#define COLS 2

PRIVATE void materialized_join(BShareTable *input1, BShareTable *input2,
                        int leftcol, int rightcol, BShareTable* result);

/**
 * Evaluates the performance of materialized join-group-by-count.
 **/

int main(int argc, char** argv) {

  if (argc < 4) {
    printf("\n\nUsage: %s <NUM_ROWS_1> <NUM_ROWS_2> <BATCH_SIZE>\n\n", argv[0]);
    return -1;
  }

  // initialize communication
  init(argc, argv);

  const int ROWS1 = atol(argv[1]); // input1 size
  const int ROWS2 = atol(argv[2]); // input2 size
  const int BATCH_SIZE_ = atol(argv[3]); // batch size

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS1, 2*COLS, 1};
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {-1, rank, ROWS2, 2*COLS, 2};
  allocate_bool_shares_table(&t2);

  // random bits required by group-by
  int prod_size = ROWS1*ROWS2;
  int num_rands = prod_size*log2(prod_size) + 1;
  AShare *ra = (AShare *) malloc(num_rands*sizeof(AShare));
  BShare *rb = (BShare *) malloc(num_rands*sizeof(BShare));

  // initialize rand bits (all equal to 1)
  for (int i=0; i<num_rands; i++) {
    ra[i] = (unsigned int) 1;
    rb[i] = (unsigned int) 1;
  }

  if (rank == 0) { //P1
    // Initialize input data and shares
    Table r1, r2;
    generate_random_table(&r1, ROWS1, COLS);
    generate_random_table(&r2, ROWS2, COLS);

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS1, 2*COLS, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS1, 2*COLS, 1};
    allocate_bool_shares_table(&t13);

    // t2 Bshare tables for P2, P3 (local to P1)
    BShareTable t22 = {-1, 1, ROWS2, 2*COLS, 2};
    allocate_bool_shares_table(&t22);
    BShareTable t23 = {-1, 2, ROWS2, 2*COLS, 2};
    allocate_bool_shares_table(&t23);
    
    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r1, &t1, &t12, &t13);
    // Generate boolean shares for r2
    generate_bool_share_tables(&r2, &t2, &t22, &t23);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS1 * 2 * COLS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t22.content[0][0]), ROWS2 * 2 * COLS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS1 * 2 * COLS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t23.content[0][0]), ROWS2 * 2 * COLS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1.content);
    free(t12.content);
    free(t13.content);
    free(r2.content);
    free(t22.content);
    free(t23.content);

  }
  else if (rank == 1) { //P2
    MPI_Recv(&(t1.content[0][0]), ROWS1 * 2 * COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS2 * 2 * COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&(t1.content[0][0]), ROWS1 * 2 * COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS2 * 2 * COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);
   
  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  /* =======================================================
     Measure naive join group-by count
  ======================================================== */
  /**
   * SELECT a, COUNT(a)
   * FROM t1, t2
   * WHERE t1.a=t2.a
   * GROUP BY a
  **/
  // start timer
  gettimeofday(&begin, 0);

  // allocate join result table
  BShareTable res_table = {-1, rank, ROWS1*ROWS2, 2*COLS+2, 1};
  allocate_bool_shares_table(&res_table);

  // Apply join first and materialize result
  materialized_join(&t1, &t2, 0, 0, &res_table);

  // we don't need original tables anymore
  free(t1.content); free(t2.content);

  // sort on group-by predicate
  unsigned int sort_att[1] = {0};
  bool asc[1] = {1};
  bitonic_sort_batch(&res_table, sort_att, 1, asc, res_table.numRows/2);

  // copy the join bits into an array to provide as an argument for
  // group-by-count
  BShare *join_selected = (BShare *) malloc(ROWS1*ROWS2*sizeof(BShare));
  assert(join_selected !=NULL);

  // convert selected bits to arithmetic
  AShare *join_selected_a = (AShare *) malloc(ROWS1*ROWS2*sizeof(AShare));
  assert(join_selected_a !=NULL);

  convert_single_bit_array(join_selected, ra, rb, ROWS1*ROWS2,
                            join_selected_a);

  // apply group-by-count on join output
  // the results are in join_selected_a
  unsigned key_index[1] = {0};
  group_by_count_sel_odd_even(&res_table, key_index, 1, BATCH_SIZE_, join_selected, join_selected_a, rb, ra);
  
  free(join_selected);
  
  // open result
  Data *open_res = (Data *) malloc(ROWS1*ROWS2*sizeof(Data));
  assert(open_res !=NULL);
  open_array(join_selected_a, ROWS1*ROWS2, open_res);
  
  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("%d\tNAIVE GROUP-BY-JOIN\t%.3f\n", ROWS1, elapsed);
  }

  free(res_table.content); free(open_res);

  // tear down communication
  MPI_Finalize();
  return 0;
}

// The result is stored in a new BShareTable whose first columns contain
// the matching pairs of the original tables and
// the last 2 columns contain the join result bits.
PRIVATE void materialized_join(BShareTable *input1, BShareTable *input2,
                        int leftcol, int rightcol, BShareTable* result) {

  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int res_index = result->numCols-2;
  BShare *temp_local = (BShare *) malloc((result->numRows)*sizeof(BShare));
  BShare *temp_remote = (BShare *) malloc((result->numRows)*sizeof(BShare));

  // compute bitwise x^y^1
  for (int i=0; i<input1->numRows; i++) {
    // copy outer input's join attribute to result table
    result->content[i][0] = input1->content[i][leftcol];
    result->content[i][1] = input1->content[i][leftcol + 1];
    for (int j=0; j<input2->numRows; j++) {
      // initialize equality
      result->content[i][res_index] = input1->content[i][leftcol] ^ input2->content[j][rightcol] ^ (~(BShare)0); // local share;
      result->content[i][res_index + 1] = input1->content[i][leftcol + 1] ^ input2->content[j][rightcol + 1] ^ (~(BShare)0); // remote share
    }
  }

  // The result is stored in the (numbits/2) rightmost bits of result, res2 elements
  for (int l=0; l<numlevels; l++) {
    for (int i=0; i<result->numRows; i++) {
      result->content[i][res_index] = eq_b_level2(numbits >> l,
                                                  result->content[i][res_index],
                                                  result->content[i][res_index + 1]);
    }

    // Exchange results of logical and, except for the final round
    // copy result column to temp_local and exchange it
    // last exchange not required by this query
    if (l < numlevels-1) {
      for (int i=0; i<result->numRows; i++) {
        temp_local[i] = result->content[i][res_index];
      }
      exchange_shares_array(temp_local, temp_remote, result->numRows);
        // copy exchanged result back to remote column
      for (int i=0; i<result->numRows; i++) {
        result->content[i][res_index + 1] = temp_remote[i];
      }
    }
  }
  free(temp_local); free(temp_remote);
}