#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define SHARE_TAG 193
#define PRIVATE static
#define COLS 2

/**
 * Evaluates the performance of optimized join-distinct
 **/

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("\n\nUsage: %s <NUM_ROWS_1> <NUM_ROWS_2>\n\n", argv[0]);
    return -1;
  }

  // initialize communication
  init(argc, argv);

  const int ROWS1 = atol(argv[1]); // input1 size
  const int ROWS2 = atol(argv[2]); // input2 size

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS1, 2*COLS, 1};
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {-1, rank, ROWS2, 2*COLS, 2};
  allocate_bool_shares_table(&t2);

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

  BShare *join_res = (BShare *) malloc(ROWS1*sizeof(BShare));
  BShare *join_remote = (BShare *) malloc(ROWS1*sizeof(BShare));

  /* =======================================================
    Measure optimized join-distinct on the same attribute
  ======================================================== */
  /**
   * SELECT DISTINCT(t1.0)
   * FROM t1, t2
   * WHERE t1.0=t2.0
  **/
  // start timer
  gettimeofday(&begin, 0);

  // sort on distinct predicate
  unsigned int sort_att[1] = {0};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, sort_att, 1, asc, t1.numRows/2);

  // apply distinct on t1
  BitShare *res_distinct = (BitShare *) malloc(ROWS1*sizeof(BitShare));
  distinct_batch(&t1, sort_att, 1, res_distinct, t1.numRows - 1);

  // exchange distinct result
  BitShare *res_distinct_remote = (BitShare *) malloc(ROWS1*sizeof(BitShare));
  exchange_bit_shares_array(res_distinct, res_distinct_remote, ROWS1);

  // join on result ==> we only need a semi-join
  in(&t1, &t2, 0, 0, join_res, t1.numRows);

  // exchange semi-join results
  exchange_shares_array(join_res, join_remote, ROWS1);

  // Open records where b_open = 1
  // b_open = b_distinct and b_join
  BShare max=0xFFFFFFFFFFFFFFFF;
  BShare *b_open = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(b_open !=NULL);
  BShare *b_open_remote = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(b_open_remote !=NULL);
  BShare *res = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(res !=NULL);

  for (int i=0; i< ROWS1; i++) { // we only need to scan the outer relation
    b_open[i] = and_b(res_distinct[i], res_distinct_remote[i],
                          join_res[i], join_remote[i],
                          get_next_rb()) & 1;
  }
  
  exchange_shares_array(b_open, b_open_remote, ROWS1);

  for (int i=0; i< ROWS1; i++) {
    BShare b1 = -b_open[i];
    BShare b2 = -b_open_remote[i];
    // res = b_open * att + (1-b) * dummy
    res[i] = and_b(b1, b2,
                   t1.content[i][0], t1.content[i][1],
                   get_next_rb());
    res[i] ^= and_b(~b1, ~b2,
                      max, max,
                      get_next_rb());
  }

  free(b_open); free(b_open_remote);

  Data *open_res = (Data *) malloc(ROWS1*sizeof(Data));
  assert(open_res !=NULL);
  
  open_b_array(res, ROWS1, open_res);
  
  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("%d\tOPT JOIN-DISTINCT\t%.3f\n", ROWS1, elapsed);
  }

  free(join_res); free(join_remote); free(res); free(open_res);
  free(res_distinct); free(res_distinct_remote); free(t1.content);

  // tear down communication
  MPI_Finalize();
  return 0;
}