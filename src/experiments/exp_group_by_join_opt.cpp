#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define SHARE_TAG 193
#define PRIVATE static
#define COLS1 6
#define COLS2 4

/**
 * Evaluates the performance of optimized join-group-by-count.
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
  const int BATCH_SIZE_ = atoi(argv[3]); // batch size for left input

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS1, 2*COLS1, 1};
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {-1, rank, ROWS2, 2*COLS2, 2};
  allocate_bool_shares_table(&t2);

  int num_rands = ROWS1*log2(ROWS1) + 1;
  BShare *rb_left = (BShare *) malloc(num_rands*sizeof(BShare));
  AShare *ra_left = (AShare *) malloc(num_rands*sizeof(AShare));
  BShare *rb_right = (BShare *) malloc(BATCH_SIZE_*ROWS2*sizeof(BShare));
  AShare *ra_right = (AShare *) malloc(BATCH_SIZE_*ROWS2*sizeof(AShare));

  // initialize rand bits (all equal to 1)
  for (int i=0; i<num_rands; i++) {
    rb_left[i] = (unsigned int) 1;
    ra_left[i] = (unsigned int) 1;
  }

  for (int i=0; i<ROWS2; i++) {
    rb_right[i] = (unsigned int) 1;
    ra_right[i] = (unsigned int) 1;
  }

  if (rank == 0) { //P1
    // Initialize input data and shares
    Table r1, r2;
    generate_random_table(&r1, ROWS1, COLS1);
    generate_random_table(&r2, ROWS2, COLS2);

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS1, 2*COLS1, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS1, 2*COLS1, 1};
    allocate_bool_shares_table(&t13);

    // t2 Bshare tables for P2, P3 (local to P1)
    BShareTable t22 = {-1, 1, ROWS2, 2*COLS2, 2};
    allocate_bool_shares_table(&t22);
    BShareTable t23 = {-1, 2, ROWS2, 2*COLS2, 2};
    allocate_bool_shares_table(&t23);
    
    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r1, &t1, &t12, &t13);
    // Generate boolean shares for r2
    generate_bool_share_tables(&r2, &t2, &t22, &t23);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t22.content[0][0]), ROWS2 * 2 * COLS2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t23.content[0][0]), ROWS2 * 2 * COLS2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1.content);
    free(t12.content);
    free(t13.content);
    free(r2.content);
    free(t22.content);
    free(t23.content);

  }
  else if (rank == 1) { //P2
    MPI_Recv(&(t1.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS2 * 2 * COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&(t1.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS2 * 2 * COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);
   
  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  /* =======================================================
     Measure optimized join group-by count
  ======================================================== */
  /**
   * SELECT a, COUNT(a)
   * FROM t1, t2
   * WHERE t1.a=t2.a
   * GROUP BY a
  **/
  // start timer
  gettimeofday(&begin, 0);

  // sort left on group-by attribute
  unsigned int sort_att[1] = {0};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, sort_att, 1, asc, t1.numRows/2);

  // apply group-by-join in batches
  int num_batches = ROWS1 / BATCH_SIZE_;
  for (int i=0; i<num_batches; i+=BATCH_SIZE_) {
    // apply group-by-join
    /*group_by_join(&t1, &t2, i, i+BATCH_SIZE_, 0, 0, 0, 2,
                rb_left, ra_left, rb_right, ra_right, 2, 4);*/
    group_by_join_first(&t1, &t2, i, i+BATCH_SIZE_, 0, 0, 2,
                        rb_right, ra_right, 2);
  }
  // Apply second phase
  unsigned key_index[1] = {0};
  group_by_sum_odd_even(&t1, BATCH_SIZE_, rb_left, ra_left, 2, 4, key_index, 1);

  Data *open_res = (Data *) malloc(ROWS1*sizeof(Data));
  assert(open_res !=NULL);
  AShare *out = (AShare *) malloc(ROWS1*sizeof(AShare));
  assert(out !=NULL);

  // reveal the group_by attribute and the aggregation (sum)
  for (int i=0; i<ROWS1; i++) {
    out[i] = t1.content[i][2];
  }

  free(t1.content);
  open_array(out, ROWS1, open_res);
  
  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("%d\tOPT group-by-join\t%.3f\n", ROWS1, elapsed);
  }

  free(out); free(open_res);
  free(ra_right); free(ra_left); free(rb_left); free(rb_right);  

  // tear down communication
  MPI_Finalize();
  return 0;
}