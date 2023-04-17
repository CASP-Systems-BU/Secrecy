#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS_C 3
#define COLS_O 4
#define C 42

/**
 * Tests the correctness of TPC-H Q13 (baseline).
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

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // 0: C_CUSTOMERKEY, 2: C_COUNT (sum), 4: count (custdist)
  BShareTable t1 = {-1, rank, ROWS_C, 2*COLS_C, 1};
  allocate_bool_shares_table(&t1);
  // 0: O_ORDERKEY, 2: O_CUSTOMERKEY, 4: O_COMMENT
  // 6: SELECTED
  BShareTable t2 = {-1, rank, ROWS_O, 2*COLS_O, 2};
  allocate_bool_shares_table(&t2);

  // shares of constant
  BShare s_comm1, s_comm2;

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

    Data orders[ROWS_O][COLS_O] = {{4, 2, 99, 0},
                                    {6, 3, 99, 0},
                                    {5, 2, 99, 0},
                                    {4, 4, 99, 0},
                                    {2, 1, 99, 0},
                                    {1, 1, 99, 0},
                                    {3, 1, 99, 0},
                                    {7, 3, 99, 0},
                                    {8, 4, 99, 0},
                                    {10, 4, C, 0},
                                    {9, 4, C, 0},
                                    {13, 4, C, 0},
                                    {15, 4, C, 0},
                                    {14, 4, 99, 0},
                                    {18, 4, 99, 0},
                                    {17, 4, 99, 0}};

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

    BShare s_comm3;

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r_customer, &t1, &t12, &t13);
    // Generate boolean shares for r2
    generate_bool_share_tables(&r_order, &t2, &t22, &t23);
    // generate shares for constant
    generate_bool_share(C, &s_comm1, &s_comm2, &s_comm3);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS_C * 2 * COLS_C, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t22.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&s_comm2, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&s_comm3, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS_C * 2 * COLS_C, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t23.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&s_comm3, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&s_comm1, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

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
    MPI_Recv(&s_comm1, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&s_comm2, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  // STEP 1: SELECTION ON ORDERS
  #if DEBUG
    if (rank==0) {
      printf("1st SELECTION.\n");
    }
  #endif

  BShare *comm = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(comm!=NULL);
  BShare *rem_comm = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_comm!=NULL);

  // populate vector
  for(int i=0; i<ROWS_O; i++) {
    comm[i] = t2.content[i][4];
    rem_comm[i] = t2.content[i][5];
  }

  BShare *sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(sel!=NULL);
  BShare *rem_sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_sel!=NULL);

  eq_b_array_const(comm, rem_comm, s_comm1, s_comm2, ROWS_O, sel);

  // compute not_b selected
  for (int i=0; i<ROWS_O; i++) {
      sel[i] ^= 1;
  }
  exchange_shares_array(sel, rem_sel, ROWS_O);

  // Copy arithmetic selection bits to columns 6, 7
  for (int i=0; i<ROWS_O; i++) {
    t2.content[i][6] = sel[i];
    t2.content[i][7] = rem_sel[i];
  }
  free(sel); free(rem_sel); free(comm); free(rem_comm);

  // STEP 2: Join CUSTOMER and ORDERS on CUSTOMERKEY (0, 2)
  BShare *join_res = (BShare *) malloc(ROWS_C*ROWS_O*sizeof(BShare));
  assert(join_res!=NULL);
  BShare *rem_join_res = (BShare *) malloc(ROWS_C*ROWS_O*sizeof(BShare));
  assert(rem_join_res!=NULL);

  Predicate_B p_join = {EQ, NULL, NULL, 0, 2};
  join_b_batch(&t1, &t2, 0, ROWS_C, 0, ROWS_O, p_join, rem_join_res, join_res);
  exchange_shares_array(join_res, rem_join_res, ROWS_C*ROWS_O);

  // STEP 3: AND join result with selection result
  for (int i=0; i<ROWS_C*ROWS_O; i++) {
    join_res[i] = and_b(join_res[i], rem_join_res[i],
                        t2.content[i % ROWS_O][6], t2.content[i % ROWS_O][7],
                        get_next_rb());
  }
  exchange_shares_array(join_res, rem_join_res, ROWS_C*ROWS_O);

  // allocate join result table
  // 0: C_CUSTOMERKEY, 2: O_CUSTOMERKEY, 4: SELECTED (c_count)
  BShareTable res_table = {-1, rank, ROWS_C*ROWS_O, 2*3, 1};
  allocate_bool_shares_table(&res_table);
  for (int i=0; i<ROWS_C; i++) {
    for (int j=0; j<ROWS_O; j++) {
      res_table.content[j + i * ROWS_O][0] = t1.content[i][0];
      res_table.content[j + i * ROWS_O][1] = t1.content[i][1];
      res_table.content[j + i * ROWS_O][2] = t2.content[j][2];
      res_table.content[j + i * ROWS_O][3] = t2.content[j][3];
      res_table.content[j + i * ROWS_O][4] = join_res[j + i * ROWS_O];
      res_table.content[j + i * ROWS_O][5] = rem_join_res[j + i * ROWS_O];
    }
  }
  
  free(join_res); free(rem_join_res);
  // we don't need original tables anymore
  free(t1.content); free(t2.content);

  // STEP 4: SORT on C_CUSTOMERKEY (0)
  #if DEBUG
    if (rank==0) {
      printf("ORDER-BY C_CUSTOMERKEY.\n");
    }
  #endif
  unsigned int att_index[1] = {0};
  bool asc[1] = {1};
  long res_len = ROWS_C*ROWS_O;
  bitonic_sort_batch(&res_table, att_index, 1, asc, res_len/2);

  // STEP 5: GROUP-BY-COUNT on C_CUSTOMERKEY (res_table.0)
  #if DEBUG
    if (rank==0) {
      printf("GROUP-BY on C_CUSTOMERKEY.\n");
    }
  #endif

  unsigned key_indices[1] = {0};
  group_by_sum_rca_odd_even(&res_table, 3, key_indices, 1);

  // STEP 6: SORT on C_COUNT (4)
  #if DEBUG
    if (rank==0) {
      printf("ORDER-BY C_COUNT.\n");
    }
  #endif
  att_index[0] = 4;
  bitonic_sort_batch(&res_table, att_index, 1, asc, res_len/2);

  /** Copy C_COUNT to cols 2, 3
   * and initialize last 2 columns to bool shares of 1 for the custdist **/
  for (int i=0; i<res_len; i++) {
    res_table.content[i][2] = res_table.content[i][4];
    res_table.content[i][3] = res_table.content[i][5];
    res_table.content[i][4] = rank % 2;
    res_table.content[i][5] = succ % 2;
  }

  // STEP 7: GROUP-BY-COUNT on C_COUNT (res_table.2)
  #if DEBUG
    if (rank==0) {
      printf("GROUP-BY on C_COUNT.\n");
    }
  #endif

  key_indices[0] = 2;
  group_by_sum_rca_odd_even(&res_table, 3, key_indices, 1);

  // STEP 8: FINAL ORDER-BY custdist desc, c_count desc
  #if DEBUG
    if (rank==0) {
      printf("ORDER-BY CUSTDIST, C_COUNT.\n");
    }
  #endif
  unsigned int ind[2] = {4, 2};
  bool asc2[2] = {0, 0};
  bitonic_sort_batch(&res_table, ind, 2, asc2, res_len/2);

  // There should be 4 valid lines in the result
  Data result[res_len][2];
  // Open all elements
  for (int i=0; i<res_len; i++) {
    result[i][0] = open_b(res_table.content[i][2]); // C_COUNT
    result[i][1] = open_b(res_table.content[i][4]); // CUSTDIST
  }

  if (rank == 0) {
    assert(result[124][0]==0);
    assert(result[124][1]==4);
    assert(result[125][0]==2);
    assert(result[125][1]==2);
    assert(result[126][0]==5);
    assert(result[126][1]==1);
    assert(result[127][0]==3);
    assert(result[127][1]==1);
    #if DEBUG
      for (int i=0; i<res_len; i++) {
        printf("[%d] (C_COUNT, CUSTDIST) = %lld, %lld\n", i, result[i][0], result[i][1]);
      }
    #endif
  }

  if (rank == 0) {
    printf("TEST TPC-H Q13 (baseline): OK.\n");
  }

  free(res_table.content);

  // tear down communication
  MPI_Finalize();
  return 0;
}
