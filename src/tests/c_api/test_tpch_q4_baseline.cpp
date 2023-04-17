#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define COLS_O 4
#define COLS_L 4
#define D1 1993
#define D2 1995

/**
 * Tests the correctness of TPC-H Q4 (baseline).
 *
 * SELECT O_ORDERPRIORITY, COUNT(*) AS ORDER_COUNT FROM ORDERS
 * WHERE O_ORDERDATE >= '1993-07-01'
 * AND O_ORDERDATE < dateadd(mm,3, cast('1993-07-01' as date))
 * AND EXISTS (
 *    SELECT * FROM LINEITEM WHERE L_ORDERKEY = O_ORDERKEY AND L_COMMITDATE < L_RECEIPTDATE)
 * GROUP BY O_ORDERPRIORITY
 * ORDER BY O_ORDERPRIORITY
 **/

int main(int argc, char** argv) {

  const long ROWS_O = 4; // ORDERS input size
  const long ROWS_L = 16; // LINEITEM input size

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // 0: O_ORDERKEY, 2: O_ORDERDATE, 4: O_ORDERPRIORITY
  // Cols 6 and 7 store the result of the selection
  BShareTable t1 = {-1, rank, ROWS_O, 2*COLS_O, 1};
  allocate_bool_shares_table(&t1);
  // 0: L_ORDERKEY, 2: L_COMMITDATE, 4: L_RECEIPTDATE
  // Cols 6 and 7 store the result of the selection
  BShareTable t2 = {-1, rank, ROWS_L, 2*COLS_L, 2};
  allocate_bool_shares_table(&t2);

  // shares of constants (dates)
  BShare sd1, sd2;

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data orders[ROWS_O][COLS_O] = {{2, 1994, 200, 0},
                                  {1, 1993, 100, 0},
                                  {4, 1994, 100, 0},
                                  {7, 1991, 300, 0}};

    Data lineitem[ROWS_L][COLS_L] = {{1, 1991, 1993, 0},
                                    {2, 1992, 1994, 0},
                                    {3, 1993, 1993, 0},
                                    {4, 1992, 1994, 0},
                                    {5, 1995, 1994, 0},
                                    {6, 1997, 1998, 0},
                                    {1, 1995, 1994, 0},
                                    {2, 1995, 1994, 0},
                                    {3, 1990, 1997, 0},
                                    {4, 2000, 1999, 0},
                                    {5, 1998, 1994, 0},
                                    {6, 1996, 1993, 0},
                                    {6, 1993, 1994, 0},
                                    {6, 1995, 1992, 0},
                                    {5, 1999, 1995, 0},
                                    {3, 1993, 1994, 0}};

    Data ** c1 = allocate_2D_data_table(ROWS_O, COLS_O);
    for (int i=0;i<ROWS_O;i++){
      for(int j=0;j<COLS_O;j++){
        c1[i][j] = orders[i][j];
      }
    }

    Data ** c2 = allocate_2D_data_table(ROWS_L, COLS_L);
    for (int i=0;i<ROWS_L;i++){
      for(int j=0;j<COLS_L;j++){
        c2[i][j] = lineitem[i][j];
      }
    }

    Table r_orders = {-1, ROWS_O, COLS_O, c1};
    Table r_lineitem = {-1, ROWS_L, COLS_L, c2};

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS_O, 2*COLS_O, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS_O, 2*COLS_O, 1};
    allocate_bool_shares_table(&t13);

    // t2 Bshare tables for P2, P3 (local to P1)
    BShareTable t22 = {-1, 1, ROWS_L, 2*COLS_L, 2};
    allocate_bool_shares_table(&t22);
    BShareTable t23 = {-1, 2, ROWS_L, 2*COLS_L, 2};
    allocate_bool_shares_table(&t23);

    BShare d12, d22, d13, d23;

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r_orders, &t1, &t12, &t13);
    // Generate boolean shares for r2
    generate_bool_share_tables(&r_lineitem, &t2, &t22, &t23);
    // generate shares for constants
    generate_bool_share(D1, &sd1, &d12, &d13);
    generate_bool_share(D2, &sd2, &d22, &d23);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t22.content[0][0]), ROWS_L * 2 * COLS_L, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&d12, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&d22, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(t23.content[0][0]), ROWS_L * 2 * COLS_L, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&d13, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&d23, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(c1);
    free(t12.content);
    free(t13.content);
    free(c2);
    free(t22.content);
    free(t23.content);

  }
  else { //P2 or_a P3
    MPI_Recv(&(t1.content[0][0]), ROWS_O * 2 * COLS_O, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(t2.content[0][0]), ROWS_L * 2 * COLS_L, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&sd1, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&sd2, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  // STEP 1: Selection L_COMMITDATE < L_RECEIPTDATE
  #if DEBUG
    if (rank==0) {
      printf("1st selection on LINEITEM.\n");
    }
  #endif

  BShare *cdate = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(cdate!=NULL);
  BShare *rem_cdate = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(rem_cdate!=NULL);

  BShare *rdate = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(rdate!=NULL);
  BShare *rem_rdate = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(rem_rdate!=NULL);

  // populate vectors
  for(int i=0; i<ROWS_L; i++) {
    cdate[i] = t2.content[i][2];
    rem_cdate[i] = t2.content[i][3];
    rdate[i] = t2.content[i][4];
    rem_rdate[i] = t2.content[i][5];
  }

  BitShare *sel_l = (BitShare *) malloc(ROWS_L*sizeof(BitShare));
  assert(sel_l!=NULL);
  BitShare *rem_sel_l = (BitShare *) malloc(ROWS_L*sizeof(BitShare));
  assert(rem_sel_l!=NULL);

  greater_batch(rdate, rem_rdate, cdate, rem_cdate, ROWS_L, sel_l);
  exchange_bit_shares_array(sel_l, rem_sel_l, ROWS_L);

  // Copy selection bits to columns 6, 7
  for (int i=0; i<ROWS_L; i++) {
    t2.content[i][6] = (BShare)sel_l[i];
    t2.content[i][7] = (BShare)rem_sel_l[i];
  }

  free(cdate); free(rem_cdate); free(rdate); free(rem_rdate);
  free(sel_l); free(rem_sel_l);

 // STEP 2: Apply selection O_ORDERDATE >= D1
  #if DEBUG
    if (rank==0) {
      printf("1st selection on ORDERS.\n");
    }
  #endif

  BShare *odate = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(odate!=NULL);
  BShare *rem_odate = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_odate!=NULL);

  // populate vector
  for(int i=0; i<ROWS_O; i++) {
    odate[i] = t1.content[i][2];
    rem_odate[i] = t1.content[i][3];
  }

  BitShare *sel = (BitShare *) malloc(ROWS_O*sizeof(BitShare));
  assert(sel!=NULL);
  BitShare *rem_sel = (BitShare *) malloc(ROWS_O*sizeof(BitShare));
  assert(rem_sel!=NULL);

  BShare rem_sd1 = exchange_shares(sd1);
  geq_batch_const(odate, rem_odate, sd1, rem_sd1, ROWS_O, sel);
  exchange_bit_shares_array(sel, rem_sel, ROWS_O);

  // STEP 3: Apply selection O_ORDERDATE < D2
  #if DEBUG
    if (rank==0) {
      printf("2nd selection on ORDERS.\n");
    }
  #endif

  BitShare *sel2 = (BitShare *) malloc(ROWS_O*sizeof(BitShare));
  assert(sel2!=NULL);
  BitShare *rem_sel2 = (BitShare *) malloc(ROWS_O*sizeof(BitShare));
  assert(rem_sel2!=NULL);

  BShare rem_sd2 = exchange_shares(sd2);
  geq_batch_const(odate, rem_odate, sd2, rem_sd2, ROWS_O, sel2);
  for(int i=0; i<ROWS_O; i++) {
    sel2[i] ^= (BitShare)1;
  }
  exchange_bit_shares_array(sel2, rem_sel2, ROWS_O);
  free(odate); free(rem_odate);

  // STEP 4: Compute AND of selections
  BShare mask = 1;
  #if DEBUG
    if (rank==0) {
      printf("AND of selections on ORDERS.\n");
    }
  #endif

  BShare *res_sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(res_sel!=NULL);
  BShare *rem_res_sel = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_res_sel!=NULL);

  for (int j=0; j<ROWS_O; j++) {
    res_sel[j] = and_b((BShare)sel[j], (BShare)rem_sel[j],
                    (BShare)sel2[j], (BShare)rem_sel2[j],
                    get_next_rb()) & mask;
  }
  exchange_shares_array(res_sel, rem_res_sel, ROWS_O);
  free(sel); free(rem_sel); free(sel2); free(rem_sel2);

  // STEP 5: Fused Selection-IN
  #if DEBUG
    if (rank==0) {
      printf("(fused) IN.\n");
    }
  #endif

  BShare *in_res = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(in_res!=NULL);
  BShare *rem_in_res = (BShare *) malloc(ROWS_O*sizeof(BShare));
  assert(rem_in_res!=NULL);
  in_sel_right(&t1, &t2, 0, 0, 6, in_res, ROWS_O);
  exchange_shares_array(in_res, rem_in_res, ROWS_O);

  #if DEBUG
    if (rank==0) {
      printf("AND with result of IN.\n");
    }
  #endif

  for (int j=0; j<ROWS_O; j++) {
    res_sel[j] = and_b(in_res[j], rem_in_res[j],
                      res_sel[j], rem_res_sel[j], get_next_rb())
                      & mask;
  }
  exchange_shares_array(res_sel, rem_res_sel, ROWS_O);

  free(in_res); free(rem_in_res);

  // Copy selected bit to the ORDERS table
  for (int i=0; i<ROWS_O; i++) {
    t1.content[i][6] = res_sel[i];
    t1.content[i][7] = rem_res_sel[i];
  }
  free(res_sel); free(rem_res_sel);

  // STEP 6: Order by O_ORDERPRIORITY (att=4)
  #if DEBUG
    if (rank==0) {
      printf("Sorting.\n");
    }
  #endif
  unsigned int att_index[1] = {4};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS_O/2);

  // STEP 7: GROUP-BY-COUNT on O_ORDERPRIORITY (att=4)
  #if DEBUG
    if (rank==0) {
      printf("Group-by.\n");
    }
  #endif

  unsigned int gr_index[1] = {4};
  group_by_sum_rca_odd_even(&t1, 3, gr_index, 1);

  // There should be 2 valid lines in the result with this order:
  // (O_ORDERPRIORITY, COUNT) = 100, 2
  // (O_ORDERPRIORITY, COUNT) = 200, 1
  // The rest of the lines in the result are garbage (due to masking)
  Data result[ROWS_O][2];
  // Open all elements
  for (int i=0; i<ROWS_O; i++) {
    result[i][0] = open_b(t1.content[i][4]); // O_ORDERPRIORITY
    result[i][1] = open_b(t1.content[i][6]); // COUNT
  }

  if (rank == 0) {
    assert(result[0][0]==100);
    assert(result[0][1]==2);
    assert(result[2][0]==200);
    assert(result[2][1]==1);
    #if DEBUG
      for (int i=0; i<ROWS_O; i++) {
        printf("[%d] (O_ORDERPRIORITY, COUNT) = %lld, %lld\n",
                 i, result[i][0], result[i][1]);
      }
    #endif
    printf("TEST TPC-H Q4 Baseline: OK.\n");
  }

  free(t1.content); free(t2.content);

  // tear down communication
  MPI_Finalize();
  return 0;
}
