#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define DEBUG 0
#define SHARE_TAG 193

// Table one: boolean
// Needed COLS1:            0:SHIPDATE, 1:QUANTITY, 2:DISCOUNT, 3:EXTENDEDPRICE
// In addition to :         4:99-SHIPDATE, 5:SHIPDATE-200, 6:QUANTITY-24, 7:DISCOUNT-0.05, 8:0.07-DISCOUNT
// In addition to SEL:      9:SEL_START, 10:SEL_END, 11:DISCOUNT_MIN, 12:DISCOUNT_MAX, 13:SEL_QUANTITY
#define COLS1 14

// Table two: arithmetic
// Needed COLS1:            0:DISCOUNT, 1:EXTENDEDPRICE
#define COLS2 2

// SELECT SUM(L_EXTENDEDPRICE*L_DISCOUNT) AS REVENUE
// FROM LINEITEM
// WHERE L_SHIPDATE >= '1994-01-01'
// AND L_SHIPDATE < dateadd(yy, 1, cast('1994-01-01' as date))
// AND L_DISCOUNT BETWEEN .06 - 0.01 AND .06 + 0.01
// AND L_QUANTITY < 24

void allocate_int_shares_table_row(AShareTable *table){
  int length = sizeof(AShare *) * table->numRows +
               sizeof(AShare) * table->numRows * table->numCols;
  table->content = (AShare **)malloc(length);
  assert(table->content != NULL);
  AShare *ptr = (AShare *)(table->content + table->numRows);
  for (int i = 0; i < table->numRows; i++)
    table->content[i] = (ptr + table->numCols * i);
}

int main(int argc, char **argv){
  if (argc < 2) {
    printf("\n\nUsage: %s <NUM_ROWS_LINEITEM>\n\n", argv[0]);
    return -1;
  }

  const int ROWS_L = atol(argv[1]);

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable tb = {-1, rank, ROWS_L, 2 * COLS1, 1};
  allocate_bool_shares_table(&tb);
  AShareTable ta = {-1, rank, ROWS_L, 2 * COLS2, 1};
  allocate_int_shares_table_row(&ta);

  BShare *rb = (BShare *) calloc(ROWS_L, sizeof(BShare));
  AShare *ra = (AShare *) calloc(ROWS_L, sizeof(AShare));

  if (rank == 0){
    // P1
    // Initialize input data and shares
    // initializing COLS1: 0:SHIPDATE - 1:QUANTITY - 2:DISCOUNT - 3:EXTENDEDPRICE
    // start date is assumed to be 100
    // assume the year end at 200

    Table r1, r2;
    generate_random_table(&r1, ROWS_L, COLS1);
    generate_random_table(&r2, ROWS_L, COLS2);

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable tb2 = {-1, 1, ROWS_L, 2 * COLS1, 1};
    allocate_bool_shares_table(&tb2);
    BShareTable tb3 = {-1, 2, ROWS_L, 2 * COLS1, 1};
    allocate_bool_shares_table(&tb3);

    // t1 Ashare tables for P2, P3 (local to P1)
    AShareTable ta2 = {-1, 1, ROWS_L, 2 * COLS2, 1};
    allocate_int_shares_table_row(&ta2);
    AShareTable ta3 = {-1, 2, ROWS_L, 2 * COLS2, 1};
    allocate_int_shares_table_row(&ta3);

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r1, &tb, &tb2, &tb3);

    // Generate arthmetic shares for r2
    generate_int_share_tables(&r2, &ta, &ta2, &ta3);

    //Send shares to P2
    MPI_Send(&(tb2.content[0][0]), ROWS_L * 2 * COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(ta2.content[0][0]), ROWS_L * 2 * COLS2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(tb3.content[0][0]), ROWS_L * 2 * COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&(ta3.content[0][0]), ROWS_L * 2 * COLS2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1.content);
    free(tb2.content);
    free(tb3.content);

    free(r2.content);
    free(ta2.content);
    free(ta3.content);
  }
  else
  {
    //P2 or P3
    MPI_Recv(&(tb.content[0][0]), ROWS_L * 2 * COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&(ta.content[0][0]), ROWS_L * 2 * COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  // start timer
  gettimeofday(&begin, 0);

// STEP 1: Check Selection Conditions
#if DEBUG
  if (rank == 0){
    printf("Evaluating Selection Conditions.\n");
  }
#endif

  // Selection Conditions
  // What operators I have: =, < 0
  // 0:SHIPDATE >= 100  or  4:99-SHIPDATE < 0
  // 0:SHIPDATE <  200  or  5:SHIPDATE-200 < 0
  // 2:QUANTITY < 24    or  6:QUANTITY-24 < 0
  // 1:DISCOUNT >= 0.05 or  Not(7:DISCOUNT-0.05 < 0)
  // 1:DISCOUNT <= 0.07 or  Not(8:0.07-DISCOUNT < 0)
  Predicate_B CONDS[5] = {{GT, NULL, NULL, 8, 8},
                          {GT, NULL, NULL, 10, 10},
                          {GT, NULL, NULL, 12, 12},
                          {GT, NULL, NULL, 14, 14},
                          {GT, NULL, NULL, 16, 16}};

  BShare *sel = (BShare *) malloc(ROWS_L * sizeof(BShare));
  assert(sel != NULL);
  BShare *rem_sel = (BShare *) malloc(ROWS_L * sizeof(BShare));
  assert(rem_sel != NULL);

  int ind = 18;
  for (int k = 0; k < 5; k++){
    select_b(tb, CONDS[k], sel);
    exchange_shares_array(sel, rem_sel, ROWS_L);

    if (k < 3){
      for (int i = 0; i < ROWS_L; i++){
        tb.content[i][ind] = sel[i];
        tb.content[i][ind + 1] = rem_sel[i];
      }
    }
    else{
      for (int i = 0; i < ROWS_L; i++){
        // we have odd number of parties so each party can just
        // xor their share with 1 to do the negation
        tb.content[i][ind] = sel[i] ^ 1;
        tb.content[i][ind + 1] = rem_sel[i] ^ 1;
      }
    }

    ind += 2;
  }

  // STEP 2: Select rows
  // Needed COLS1:             0:SHIPDATE, 1:QUANTITY, 2:DISCOUNT, 3:EXTENDEDPRICE
  // In addition to :         4:99-SHIPDATE, 5:SHIPDATE-200, 6:QUANTITY-24, 7:DISCOUNT-0.05, 8:0.07-DISCOUNT
  // In addition to SEL:      9:SEL_START, 10:SEL_END, 11:DISCOUNT_MIN, 12:DISCOUNT_MAX, 13:SEL_QUANTITY
  // Selection order:
  // AND(AND(AND(1, 2), AND(3,4))), 5)
#if DEBUG
  if (rank == 0){
    printf("Applying Selection.\n");
  }
#endif

  // 1, 2 >> 2
  and_b_table(tb, 18, 20, ROWS_L, sel);
  exchange_shares_array(sel, rem_sel, ROWS_L);
  for (int i = 0; i < ROWS_L; i++){
    tb.content[i][20] = sel[i];
    tb.content[i][21] = rem_sel[i];
  }

  // 3, 4 >> 4
  and_b_table(tb, 22, 24, ROWS_L, sel);
  exchange_shares_array(sel, rem_sel, ROWS_L);
  for (int i = 0; i < ROWS_L; i++){
    tb.content[i][24] = sel[i];
    tb.content[i][25] = rem_sel[i];
  }

  // 2, 5 >> 5
  and_b_table(tb, 18, 26, ROWS_L, sel);
  exchange_shares_array(sel, rem_sel, ROWS_L);
  for (int i = 0; i < ROWS_L; i++){
    tb.content[i][26] = sel[i];
    tb.content[i][27] = rem_sel[i];
  }

  // 4, 5 >> 5
  and_b_table(tb, 24, 26, ROWS_L, sel);
  exchange_shares_array(sel, rem_sel, ROWS_L);
  for (int i = 0; i < ROWS_L; i++){
    tb.content[i][26] = sel[i];
    tb.content[i][27] = rem_sel[i];
  }

  // STEP 3: Change boolean to arth
#if DEBUG
  if (rank == 0){
    printf("Changing Selection boolean share to arithmetic.\n");
  }
#endif
  AShare *res = (AShare *) malloc(ROWS_L * sizeof(AShare));
  assert(res != NULL);
  convert_single_bit_array(sel, ra, rb, ROWS_L, res);
  free(sel);
  free(rem_sel);

  AShare *res_remote = (AShare *) malloc(ROWS_L * sizeof(AShare));
  assert(res_remote != NULL);
  exchange_shares_array(res, res_remote, ROWS_L);

  // STEP 4: Multiplication
#if DEBUG
  if (rank == 0){
    printf("Changing Selection boolean share to arithmetic.\n");
  }
#endif
  AShare *mulRes = (AShare *) malloc(ROWS_L * sizeof(AShare));
  assert(mulRes != NULL);

  AShare *mulResRemote = (AShare *) malloc(ROWS_L * sizeof(AShare));
  assert(mulResRemote != NULL);

  // 2:DISCOUNT * 3:EXTENDEDPRICE
  for (int i = 0; i < ROWS_L; i++){
    mulRes[i] = mul(ta.content[i][0], ta.content[i][1], ta.content[i][2], ta.content[i][3], get_next_r());
  }
  exchange_shares_array(mulRes, mulResRemote, ROWS_L);

  long long summation = 0;
  // choosing those with selection line of 1
  for (int i = 0; i < ROWS_L; i++){
    mulRes[i] = mul(mulRes[i], mulResRemote[i], res[i], res_remote[i], get_next_r());
    summation += mulRes[i];
  }

  // Step Four: getting arthmetic shares
  AShare summation_share = summation;
  open_a(summation_share);

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro * 1e-6;

  if(rank == 0) {
    printf("Time: %lf\n", elapsed);
  }

  free(tb.content);
  free(ta.content);

  // tear down communication
  MPI_Finalize();
  return 0;
}
