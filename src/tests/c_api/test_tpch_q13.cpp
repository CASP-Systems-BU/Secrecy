#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS_C 3
#define COLS_O 5
#define C 42

/**
 * Tests the correctness of TPC-H Q13.
 *
 * select c_count, count(*) as custdist
 * from (
 *  select c_custkey, count(o_orderkey)
 *    from customer left outer join orders 
 *    on c_custkey = o_custkey
 *    and o_comment <> 'C'
 *    group by c_custkey
 *  )as c_orders (c_custkey, c_count)
 * group by c_count
 * order by custdist desc, c_count desc;
 **/

int main(int argc, char** argv) {

  const long ROWS_C = 8; // CUSTOMER input size
  const long ROWS_O = 16; // ORDERS input size
  const int BATCH_SIZE_ = 2; // left batch size for semi-join

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // 0: C_CUSTOMERKEY, 2: C_COUNT (sum), 4: count (custdist)
  BShareTable t1 = {-1, rank, ROWS_C, 2*COLS_C, 1};
  allocate_bool_shares_table(&t1);
  // 0: O_ORDERKEY, 2: O_CUSTOMERKEY, 4: O_COMMENT-C,
  // 6: C-O_COMMENT, 8: SELECTED
  BShareTable t2 = {-1, rank, ROWS_O, 2*COLS_O, 2};
  allocate_bool_shares_table(&t2);

  // initialize rand bits for conversion (all equal to 0)
  int num_rands = ROWS_C*log2(ROWS_C) + 1;
  BShare *rb_left = (BShare *) calloc(num_rands, sizeof(BShare));
  AShare *ra_left = (AShare *) calloc(num_rands, sizeof(AShare));
  // initialize rand bits for group-by (all equal to 0)
  BShare *rb_right = (BShare *) calloc(ROWS_O*BATCH_SIZE_, sizeof(BShare));
  AShare *ra_right = (AShare *) calloc(ROWS_O*BATCH_SIZE_, sizeof(AShare));

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data customer[ROWS_C][COLS_C] = {{2, 0, 0},
                                  {1, 0, 0},
                                  {3, 0, 0},
                                  {5, 0, 0},
                                  {6, 0, 0},
                                  {8, 0, 0},
                                  {7, 0, 0},
                                  {4, 0, 0}};

    Data orders[ROWS_O][COLS_O] = {{4, 2, 99-C, C-99, 0},
                                    {6, 3, 99-C, C-99, 0},
                                    {5, 2, 99-C, C-99, 0},
                                    {4, 4, 99-C, C-99, 0},
                                    {2, 1, 99-C, C-99, 0},
                                    {1, 1, 99-C, C-99, 0},
                                    {3, 1, 99-C, C-99, 0},
                                    {7, 3, 99-C, C-99, 0},
                                    {8, 4, 99-C, C-99, 0},
                                    {10, 4, 0, 0, 0},
                                    {9, 4, 0, 0, 0},
                                    {13, 4, 0, 0, 0},
                                    {15, 4, 0, 0, 0},
                                    {14, 4, 99-C, C-99, 0},
                                    {18, 4, 99-C, C-99, 0},
                                    {17, 4, 99-C, C-99, 0}};

    Data ** c1 = allocate_2D_data_table(ROWS_C, COLS_C);
    for (int i=0;i<ROWS_C;i++){
      for(int j=0;j<COLS_C;j++){
        c1[i][j] = customer[i][j];
      }
    }

    Data ** c2 = allocate_2D_data_table(ROWS_O, COLS_O);
    for (int i=0;i<ROWS_O;i++){
      for(int j=0;j<COLS_O;j++){
        c2[i][j] = orders[i][j];
      }
    }

    Table r_customer = {-1, ROWS_C, COLS_C, c1};
    Table r_order = {-1, ROWS_O, COLS_O, c2};

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS_C, 2*COLS_C, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS_C, 2*COLS_C, 1};
    allocate_bool_shares_table(&t13);

    // t2 Bshare tables for P2, P3 (local to P1)
    BShareTable t22 = {-1, 1, ROWS_O, 2*COLS_O, 2};
    allocate_bool_shares_table(&t22);
    BShareTable t23 = {-1, 2, ROWS_O, 2*COLS_O, 2};
    allocate_bool_shares_table(&t23);

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r_customer, &t1, &t12, &t13);
    // Generate boolean shares for r2
    generate_bool_share_tables(&r_order, &t2, &t22, &t23);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS_C * 2 * COLS_C, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t22.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS_C * 2 * COLS_C, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t23.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(c1);
    free(t12.content);
    free(t13.content);
    free(c2);
    free(t22.content);
    free(t23.content);

  }
  else { //P2 or_a P3
    MPI_Recv(&(t1.content[0][0]), ROWS_C * 2 * COLS_C, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  // STEP 1: SORT CUSTOMER on C_CUSTOMERKEY (0)
  #if DEBUG
    if (rank==0) {
      printf("ORDER-BY C_CUSTOMERKEY.\n");
    }
  #endif
  unsigned int att_index[1] = {0};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS_C/2);

  // STEP 2: SELECTION ON ORDERS
  #if DEBUG
    if (rank==0) {
      printf("1st SELECTION.\n");
    }
  #endif

  BShare *sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(sel!=NULL);

  Predicate_B p = {EQ, NULL, NULL, 4, 6};
  select_b(t2, p, sel);

  // compute not_b selected
  for (int i=0; i<ROWS_O; i++) {
      sel[i] ^= 1;
  }

  /** Convert selected bits to arithmetic **/
  AShare *sel_a = (AShare *) malloc(ROWS_O*sizeof(AShare));
  assert(sel_a!=NULL);
  AShare *rem_sel_a = (AShare *) malloc(ROWS_O*sizeof(AShare));
  assert(rem_sel_a!=NULL);

  convert_single_bit_array(sel, ra_right, rb_right, ROWS_O, sel_a);
  exchange_a_shares_array(sel_a, rem_sel_a, ROWS_O);

  free(sel);

  // Copy arithmetic selection bits to columns 8, 9
  for (int i=0; i<ROWS_O; i++) {
    t2.content[i][8] = sel_a[i];
    t2.content[i][9] = rem_sel_a[i];
  }

  free(sel_a); free(rem_sel_a);

  // STEP 3: Fused group-by-cnt-join
  #if DEBUG
    if (rank==0) {
      printf("GROUP-BY & JOIN.\n");
    }
  #endif

  // init sum and count columns of the left table
  for (int i=0; i<ROWS_C; i++) {
    t1.content[i][2] = 0;
    t1.content[i][3] = 0;
    t1.content[i][4] = rank % 2;
    t1.content[i][5] = succ % 2;
  }

  // apply group-by-join in batches
  for (int i=0; i<ROWS_C; i+=BATCH_SIZE_) {
    // apply group-by-join
    group_by_join_first(&t1, &t2, i, i+BATCH_SIZE_, 0, 2, 8,
                        rb_right, ra_right, 2);
  }
  // apply second phase
  unsigned key_index[1] = {0};
  group_by_sum_odd_even(&t1, 3, rb_left, ra_left, 2, 4, key_index, 1);

  free(ra_left); free(rb_left); free(ra_right); free(rb_right);

  AShare *c_counts = (AShare *) malloc(ROWS_C*sizeof(AShare));
  AShare *c_remote_counts = (AShare *) malloc(ROWS_C*sizeof(AShare));
  BShare *c_counts_b = (BShare *) malloc(ROWS_C*sizeof(BShare));
  BShare *c_remote_counts_b = (BShare *) malloc(ROWS_C*sizeof(BShare));

  for (int i=0; i<ROWS_C; i++) {
    c_counts[i] = t1.content[i][2];
  }

  exchange_a_shares_array(c_counts, c_remote_counts, ROWS_C);

  // Convert c_count to binary **/
  convert_a_to_b_array(c_counts, c_remote_counts, c_counts_b, c_remote_counts_b, ROWS_C);

  // Copy binary c_counts to columns 2, 3
  for (int i=0; i<ROWS_C; i++) {
    t1.content[i][2] = c_counts_b[i];
    t1.content[i][3] = c_remote_counts_b[i];
  }

  free(c_counts); free(c_remote_counts); free(c_counts_b); free(c_remote_counts_b);

  // STEP 4: SORT ON C_COUNT (t1.2)
  #if DEBUG
    if (rank==0) {
      printf("ORDER-BY C_COUNT.\n");
    }
  #endif
  att_index[0] = 2;
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS_C/2);

  // STEP 5: GROUP-BY-COUNT on C_COUNT (t1.2)
  #if DEBUG
    if (rank==0) {
      printf("GROUP-BY on C_COUNT.\n");
    }
  #endif

  AShare *rand_a = (AShare *) calloc(num_rands, sizeof(AShare));
  BShare *rand_b = (BShare *) calloc(num_rands, sizeof(BShare));
  AShare *counters = (AShare *) malloc(ROWS_C*sizeof(AShare));
  AShare *remote_counters = (AShare *) malloc(ROWS_C*sizeof(AShare));

  // initialize counters
  for (int i=0; i<ROWS_C; i++) {
    counters[i] = rank % 2;
    remote_counters[i] = succ % 2;
  }

  unsigned key_indices[1] = {2};
  group_by_count_odd_even(&t1, key_indices, 1, 3, counters, remote_counters,
                          rand_b, rand_a);

  free(rand_a); free(rand_b);

  BShare *counters_b = (BShare *) malloc(ROWS_C*sizeof(BShare));
  BShare *remote_counters_b = (BShare *) malloc(ROWS_C*sizeof(BShare));

  // Convert custdist to binary
  convert_a_to_b_array(counters, remote_counters, counters_b, remote_counters_b, ROWS_C);

  free(counters); free(remote_counters);

  // copy custdist (counters) to columns 4, 5
  for (int i=0; i<ROWS_C; i++) {
    t1.content[i][4] = counters_b[i];
    t1.content[i][5] = remote_counters_b[i];
  }

  // STEP 6: FINAL ORDER-BY custdist desc, c_count desc
  #if DEBUG
    if (rank==0) {
      printf("ORDER-BY CUSTDIST, C_COUNT.\n");
    }
  #endif
  unsigned int ind[2] = {4, 2};
  bool asc2[2] = {0, 0};
  bitonic_sort_batch(&t1, ind, 2, asc2, ROWS_C/2);

  // There should be 4 valid lines in the result
  Data result[ROWS_C][2];
  // Open all elements
  for (int i=0; i<ROWS_C; i++) {
    result[i][0] = open_b(t1.content[i][2]); // C_COUNT
    result[i][1] = open_b(t1.content[i][4]); // CUSTDIST
  }

  if (rank == 0) {
    assert(result[4][0]==0);
    assert(result[4][1]==4);
    assert(result[5][0]==2);
    assert(result[5][1]==2);
    assert(result[6][0]==5);
    assert(result[6][1]==1);
    assert(result[7][0]==3);
    assert(result[7][1]==1);
    #if DEBUG
      for (int i=0; i<ROWS_C; i++) {
        printf("[%d] (C_COUNT, CUSTDIST) = %lld, %lld\n", i, result[i][0], result[i][1]);
      }
    #endif
  }

  if (rank == 0) {
    printf("TEST TPC-H Q13: OK.\n");
  }

  free(t1.content); free(t2.content);

  // tear down communication
  MPI_Finalize();
  return 0;
}
