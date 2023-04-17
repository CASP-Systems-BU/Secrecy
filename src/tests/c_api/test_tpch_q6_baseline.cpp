#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define COLS_L 5
#define D1 1994
#define D2 1996
#define DISC_1 1
#define DISC_2 6
#define QUANT 24

/**
 * Tests the correctness of TPC-H Q6 (baseline).
 *
 * SELECT SUM(L_EXTENDEDPRICE*L_DISCOUNT) AS REVENUE
 * FROM LINEITEM
 * WHERE L_SHIPDATE >= '1994-01-01' 
 * AND L_SHIPDATE < dateadd(yy, 1, cast('1994-01-01' as date))
 * AND L_DISCOUNT BETWEEN .06 - 0.01 AND .06 + 0.01 
 * AND L_QUANTITY < 24
 **/

int main(int argc, char** argv) {

  const long ROWS_L = 16; // LINEITEM input size

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // 0: L_EXTENDEDPRICE, 2: L_DISCOUNT, 4: L_SHIPDATE,
  // 6: L_QUANTITY
  // Cols 8 and 9 store the result of the selection
  AShareTable t1 = {-1, rank, ROWS_L, 2*COLS_L, 2};
  allocate_a_shares_table(&t1);

  // shares of constants (dates)
  AShare sd1, sd2, disc1, disc2, q;

  if (rank == 0) { //P1
    // Initialize input data and shares

    Data lineitem[ROWS_L][COLS_L] = {{1, 2, 1991, 23, 0},
                                    {2, 0, 1992, 23, 0},
                                    {3, 0, 1993, 23, 0},
                                    {4, 0, 1992, 24, 0},
                                    {5, 2, 1995, 22, 0},
                                    {6, 0, 1997, 24, 0},
                                    {1, 0, 1995, 24, 0},
                                    {2, 0, 1995, 25, 0},
                                    {3, 0, 1990, 25, 0},
                                    {4, 0, 2000, 23, 0},
                                    {5, 2, 1998, 23, 0},
                                    {6, 2, 1996, 23, 0},
                                    {6, 2, 1993, 24, 0},
                                    {6, 0, 1995, 24, 0},
                                    {5, 0, 1999, 25, 0},
                                    {3, 0, 1993, 25, 0}};

    Data ** c = allocate_2D_data_table(ROWS_L, COLS_L);
    for (int i=0;i<ROWS_L;i++){
      for(int j=0;j<COLS_L;j++){
        c[i][j] = lineitem[i][j];
      }
    }

    Table r_lineitem = {-1, ROWS_L, COLS_L, c};

    // t2 Bshare tables for P2, P3 (local to P1)
    AShareTable t12 = {-1, 1, ROWS_L, 2*COLS_L, 2};
    allocate_a_shares_table(&t12);
    AShareTable t13 = {-1, 2, ROWS_L, 2*COLS_L, 2};
    allocate_a_shares_table(&t13);

    AShare d12, d22, d13, d23, disc12, disc22, disc13, disc23, q2, q3;

    init_sharing();

    // Generate shares for r1
    // NOTE: we use arithmetic sharing
    generate_int_share_tables(&r_lineitem, &t1, &t12, &t13);

    // generate shares for constants
    generate_int_share(D1, &sd1, &d12, &d13);
    generate_int_share(D2, &sd2, &d22, &d23);
    generate_int_share(DISC_1, &disc1, &disc12, &disc13);
    generate_int_share(DISC_2, &disc2, &disc22, &disc23);
    generate_int_share(QUANT, &q, &q2, &q3);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS_L * 2 * COLS_L, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&d12, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&d22, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&disc12, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&disc22, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&q2, 1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS_L * 2 * COLS_L, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&d13, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&d23, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&disc13, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&disc23, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&q3, 1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(c);
    free(t12.content);
    free(t13.content);
  }
  else { //P2 or_a P3
    MPI_Recv(&(t1.content[0][0]), ROWS_L * 2 * COLS_L, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&sd1, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&sd2, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&disc1, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&disc2, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&q, 1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  // STEP 1: Selection L_SHIPDATE >= D1 (col 4)
  #if DEBUG
    if (rank==0) {
      printf("1st selection on LINEITEM.\n");
    }
  #endif

  BShare *a_ldate = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(a_ldate!=NULL);
  BShare *a_rem_ldate = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(a_rem_ldate!=NULL);

  BShare *lrec = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(lrec!=NULL);
  BShare *rem_lrec = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(rem_lrec!=NULL);

  BitShare *sel = (BitShare *) malloc(ROWS_L*sizeof(BShare));
  assert(sel!=NULL);
  BitShare *rem_sel = (BitShare *) malloc(ROWS_L*sizeof(BShare));
  assert(rem_sel!=NULL);

  // populate vector
  for(int i=0; i<ROWS_L; i++) {
    a_ldate[i] = t1.content[i][4];
    a_rem_ldate[i] = t1.content[i][5];
  }

  // convert L_SHIPDATE to boolean shares
  convert_a_to_b_array(a_ldate, a_rem_ldate, lrec, rem_lrec, ROWS_L);

  free(a_ldate); free(a_rem_ldate);

  /** convert constant to boolean **/
  AShare rem_sd1 = exchange_shares(sd1);
  BShare sd1_b, rem_sd1_b;
  convert_a_to_b_array(&sd1, &rem_sd1, &sd1_b, &rem_sd1_b, 1);

  geq_batch_const(lrec, rem_lrec, sd1_b, rem_sd1_b, ROWS_L, sel);
  exchange_bit_shares_array(sel, rem_sel, ROWS_L);

  // Copy selection bits to columns 8, 9
  for (int i=0; i<ROWS_L; i++) {
    t1.content[i][8] = (BShare)sel[i];
    t1.content[i][9] = (BShare)rem_sel[i];
  }

 // STEP 2: Apply selection L_SHIPDATE < D2 (col 4)
 // computed as ~ (L_SHIPDATE >= D2)
  #if DEBUG
    if (rank==0) {
      printf("2nd selection on LINEITEM.\n");
    }
  #endif

  /** convert constant to boolean **/
  AShare rem_sd2 = exchange_shares(sd2);
  BShare sd2_b, rem_sd2_b;
  convert_a_to_b_array(&sd2, &rem_sd2, &sd2_b, &rem_sd2_b, 1);
  geq_batch_const(lrec, rem_lrec, sd2_b, rem_sd2_b, ROWS_L, sel);
  // compute not_b selected
  for (int i=0; i<ROWS_L; i++) {
      sel[i] ^= 1;
  }
  exchange_bit_shares_array(sel, rem_sel, ROWS_L);

  // and with previous selection
   for (int i=0; i<ROWS_L; i++) {
     sel[i] = and_b((BShare)sel[i], (BShare)rem_sel[i],
                    t1.content[i][8], t1.content[i][9], get_next_rb()) & 1;
  }
  exchange_bit_shares_array(sel, rem_sel, ROWS_L);
  
  // Copy selection bits to columns 8, 9
  for (int i=0; i<ROWS_L; i++) {
    t1.content[i][8] = (BShare)sel[i];
    t1.content[i][9] = (BShare)rem_sel[i];
  }

  // STEP 3: Apply selection L_DISCOUNT > DISC_1 (col 2)
  #if DEBUG
    if (rank==0) {
      printf("3rd selection on LINEITEM.\n");
    }
  #endif

  BShare *a_ldisc = (BShare *)malloc(ROWS_L*sizeof(BShare));
  assert(a_ldisc!=NULL);
  BShare *a_rem_ldisc = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(a_rem_ldisc!=NULL);

  // populate vector
  for(int i=0; i<ROWS_L; i++) {
    a_ldisc[i] = t1.content[i][2];
    a_rem_ldisc[i] = t1.content[i][3];
  }

  // convert L_DISC to boolean shares
  convert_a_to_b_array(a_ldisc, a_rem_ldisc, lrec, rem_lrec, ROWS_L);

  /** convert constant to boolean **/
  AShare rem_disc1 = exchange_shares(disc1);
  BShare disc1_b, rem_disc1_b;
  convert_a_to_b_array(&disc1, &rem_disc1, &disc1_b, &rem_disc1_b, 1);

  greater_batch_const(lrec, rem_lrec, disc1_b, rem_disc1_b, ROWS_L, sel);
  exchange_bit_shares_array(sel, rem_sel, ROWS_L);

  // and with previous selection
   for (int i=0; i<ROWS_L; i++) {
     sel[i] = and_b((BShare)sel[i], (BShare)rem_sel[i],
                    t1.content[i][8], t1.content[i][9], get_next_rb()) & 1;
  }
  exchange_bit_shares_array(sel, rem_sel, ROWS_L);
  
  // Copy selection bits to columns 8, 9
  for (int i=0; i<ROWS_L; i++) {
    t1.content[i][8] = (BShare)sel[i];
    t1.content[i][9] = (BShare)rem_sel[i];
  }

 // STEP 4: Apply selection L_DISCOUNT < DISC_2 (col 2)
 // computed as ~ (L_DISCOUNT >= D2)
  #if DEBUG
    if (rank==0) {
      printf("4th selection on LINEITEM.\n");
    }
  #endif

  /** convert constant to boolean **/
  AShare rem_disc2 = exchange_shares(disc2);
  BShare disc2_b, rem_disc2_b;
  convert_a_to_b_array(&disc2, &rem_disc2, &disc2_b, &rem_disc2_b, 1);
  
  geq_batch_const(lrec, rem_lrec, disc2_b, rem_disc2_b, ROWS_L, sel);
  // compute not_b selected
  for (int i=0; i<ROWS_L; i++) {
      sel[i] ^= 1;
  }
  exchange_bit_shares_array(sel, rem_sel, ROWS_L);

  // and with previous selection
   for (int i=0; i<ROWS_L; i++) {
     sel[i] = and_b((BShare)sel[i], (BShare)rem_sel[i],
                    t1.content[i][8], t1.content[i][9], get_next_rb()) & 1;
  }

  exchange_bit_shares_array(sel, rem_sel, ROWS_L);
  
  // Copy selection bits to columns 8, 9
  for (int i=0; i<ROWS_L; i++) {
    t1.content[i][8] = (BShare)sel[i];
    t1.content[i][9] = (BShare)rem_sel[i];
  }

 // STEP 5: Apply selection L_QUANTITY < Q (col 6)
 // computed as ~ (L_QUANTITY >= Q)
  #if DEBUG
    if (rank==0) {
      printf("5th selection on LINEITEM.\n");
    }
  #endif

  // populate vector
  for(int i=0; i<ROWS_L; i++) {
    a_ldisc[i] = t1.content[i][6];
    a_rem_ldisc[i] = t1.content[i][7];
  }

  // convert L_DISC to boolean shares
  convert_a_to_b_array(a_ldisc, a_rem_ldisc, lrec, rem_lrec, ROWS_L);

   /** convert constant to boolean **/
  AShare rem_q = exchange_shares(q);
  BShare q_b, rem_q_b;
  convert_a_to_b_array(&q, &rem_q, &q_b, &rem_q_b, 1);

  geq_batch_const(lrec, rem_lrec, q_b, rem_q_b, ROWS_L, sel);
  // compute not_b selected
  for (int i=0; i<ROWS_L; i++) {
      sel[i] ^= 1;
  }
  exchange_bit_shares_array(sel, rem_sel, ROWS_L);

  // and with previous selection
  for (int i=0; i<ROWS_L; i++) {
     sel[i] = and_b((BShare)sel[i], (BShare)rem_sel[i],
                    t1.content[i][8], t1.content[i][9], get_next_rb()) & 1;
  }

  exchange_bit_shares_array(sel, rem_sel, ROWS_L);
  free(a_ldisc); free(a_rem_ldisc);

  // STEP 6: Convert selected from boolean to arithmetic
  #if DEBUG
    if (rank==0) {
      printf("Conversion.\n");
    }
  #endif

  AShare *a_sel = (AShare *) malloc(ROWS_L*sizeof(AShare));
  assert(a_sel!=NULL);
  AShare *rem_a_sel = (AShare *) malloc(ROWS_L*sizeof(AShare));
  assert(rem_a_sel!=NULL);

  AShare *ra = (AShare *) calloc(ROWS_L, sizeof(AShare));
  assert(ra!=NULL);
  BShare *rb = (BShare *) calloc(ROWS_L, sizeof(BShare));
  assert(rb!=NULL);

  BShare *bshare_sel = (BShare *) malloc(ROWS_L*sizeof(BShare));
  assert(bshare_sel!=NULL);

  for (int i=0; i<ROWS_L; i++) {
      bshare_sel[i] = sel[i];
  }

  convert_single_bit_array(bshare_sel, ra, rb, ROWS_L, a_sel);
  exchange_a_shares_array(a_sel, rem_a_sel, ROWS_L);

  free(sel); free(bshare_sel); free(rem_sel);
  free(ra); free(rb);

  // STEP 7: Multiplication L_EXTENDEDPRICE*L_DISCOUNT (0, 2) and sum
  AShare final_sum = 0;
  AShare *sum = (AShare *) calloc(ROWS_L, sizeof(AShare));
  assert(sum!=NULL);
  AShare *rem_sum = (AShare *) calloc(ROWS_L, sizeof(AShare));
  assert(rem_sum!=NULL);

  for (int i=0; i<ROWS_L; i++) {
    sum[i] = mul(t1.content[i][0], t1.content[i][1],
                 t1.content[i][2], t1.content[i][3], get_next_r());
  }

  exchange_a_shares_array(sum, rem_sum, ROWS_L);

  for (int i=0; i<ROWS_L; i++) {
    final_sum += mul(a_sel[i], rem_a_sel[i],
                    sum[i], rem_sum[i], get_next_r());
  }

  free(a_sel); free(rem_a_sel);

  // Open sum
  Data result = open_a(final_sum);
  
  if (rank == 0) {
    assert(result==10);
    #if DEBUG
      printf("sum = %lld\n", result);
    #endif
    printf("TEST TPC-H Q6 Baseline: OK.\n");
  }

  free(t1.content);

  // tear down communication
  MPI_Finalize();
  return 0;
}
