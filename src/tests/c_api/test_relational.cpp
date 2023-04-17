#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define SHARE_TAG 193
#define ROWS1 6
#define COLS1 4
#define ROWS2 4
#define COLS2 4
#define EQ1_CONST 99
#define EQ2_CONST -4

/*******************************************************************
 * Test multiple relational operators.
 *
 *  Input Relations
 * -------------------
 *      r1            r2
 *  id  |   att      id |   att
 *  1       99       1      -4
 *  2       -8       3      2
 *  3       99       5      -4
 *  4       99       7      -4
 *  5       99
 *  6       14
 *
 *  QUERY #1:
 * ------------
 *  SELECT COUNT(*)
 *  FROM r1, r2 on r1.id = r2.id
 *  WHERE r1.att = 99 AND r2.att = -4
 *
 *  Query #1 result: 2
 *
 * **************************/
int main(int argc, char** argv) {

  /** ============================================================= *
   *  INITIALIZE COMMUNICATION
   *  ============================================================= */
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  /** ============================================================= *
   *  DEFINE INTERMEDIATE RESULT AND HELPER VARIABLES
   *  ============================================================= */
  BShare rb[ROWS1*ROWS2]; // random bit binary shares
  AShare ra[ROWS1*ROWS2]; // random bit arithmetic shares

  BShare s_res1[ROWS1+1], s_res2[ROWS2]; // selection results
  //remote shares for selection results
  BShare s2_res1[ROWS1+1], s2_res2[ROWS2];

  /** ============================================================= *
   *  PERFORM PRE-PROCESSING OPERATIONS
   *  ============================================================= */
  exchange_rsz_seeds(succ, pred);

  // Exchange ROWS1*ROWS2 random bit shares
  if (rank == 0) {

    init_sharing();

    BShare rb2[ROWS1*ROWS2], rb3[ROWS1*ROWS2];
    AShare ra2[ROWS1*ROWS2], ra3[ROWS1*ROWS2];

    // Generate random bits and corresponding shares
    generate_rand_bit_shares(rb, ra, rb2, ra2, rb3, ra3, ROWS1*ROWS2);

    // Send random bit shares
    // ...to P2
    MPI_Send(&rb2, ROWS1*ROWS2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ra2, ROWS1*ROWS2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    // ...to P3
    MPI_Send(&rb3, ROWS1*ROWS2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ra3, ROWS1*ROWS2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

  }
  else if (rank == 1) { //P2
    MPI_Recv(&rb, ROWS1*ROWS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ra, ROWS1*ROWS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&rb, ROWS1*ROWS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ra, ROWS1*ROWS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  /** ============================================================= *
   *  INITIALIZE INPUT DATA AND SEND SHARES
   *  ============================================================= */

  // r1: first relation, r2: second relation
  BShare r1s1[ROWS1][COLS1], r1s2[ROWS1][COLS1], r1s3[ROWS1][COLS1],
         r2s1[ROWS2][COLS2], r2s2[ROWS2][COLS2], r2s3[ROWS2][COLS2];

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data r1[ROWS1][COLS1] = {{1, 99, 99-EQ1_CONST, EQ1_CONST-99},
                            {2, -8, -8-EQ1_CONST, EQ1_CONST+8},
                            {3, 99, 99-EQ1_CONST, EQ1_CONST-99},
                            {4, 99, 99-EQ1_CONST, EQ1_CONST-99},
                            {5, 99, 99-EQ1_CONST, EQ1_CONST-99},
                            {6, 14, 14-EQ1_CONST, EQ1_CONST-14}};

    Data r2[ROWS2][COLS2] = {{1, -4, -4-EQ2_CONST, EQ2_CONST+4},
                            {3, 2, 2-EQ2_CONST, EQ2_CONST-2},
                            {5, -4, -4-EQ2_CONST, EQ2_CONST+4},
                            {7, -4, -4-EQ2_CONST, EQ2_CONST+4}};

    // generate r1 shares
    for (int i=0; i<ROWS1; i++) {
        for (int j=0; j<COLS1; j++) {
            generate_bool_share(r1[i][j], &r1s1[i][j], &r1s2[i][j], &r1s3[i][j]);
        }
    }

    // generate r2 shares
    for (int i=0; i<ROWS2; i++) {
        for (int j=0; j<COLS2; j++) {
            generate_bool_share(r2[i][j], &r2s1[i][j], &r2s2[i][j], &r2s3[i][j]);
        }
    }

    //Send shares to P2
    MPI_Send(&r1s2[0][0], ROWS1*COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s2[0][0], ROWS2*COLS2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r1s3[0][0], ROWS1*COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s3[0][0], ROWS2*COLS2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&r1s3[0][0], ROWS1*COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s3[0][0], ROWS2*COLS2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r1s1[0][0], ROWS1*COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r2s1[0][0], ROWS2*COLS2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&r1s1[0][0], ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r2s1[0][0], ROWS2*COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r1s2[0][0], ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r2s2[0][0], ROWS2*COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&r1s1[0][0], ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r2s1[0][0], ROWS2*COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r1s2[0][0], ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r2s2[0][0], ROWS2*COLS2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  // allocate BShare tables
  BShareTable in1 = {-1, rank, ROWS1, 2*COLS1, 1}, in2 = {-1, rank, ROWS2, 2*COLS2, 2};
  allocate_bool_shares_table(&in1);
  allocate_bool_shares_table(&in2);

  // copy shares into the BShareTables
  for (int i=0; i<ROWS1; i++) {
    for (int j=0; j<COLS1; j++) {
        in1.content[i][2 * j] = r1s1[i][j];
        in1.content[i][2 * j + 1] = r1s2[i][j];
    }
  }
  // relation 2
  for (int i=0; i<ROWS2; i++) {
    for (int j=0; j<COLS2; j++) {
        in2.content[i][2 * j] = r2s1[i][j];
        in2.content[i][2 * j + 1] = r2s2[i][j];
    }
  }

  /** ============================================================= *
   *  COMPUTE THE QUERY
   *  ============================================================= */

  // 1. Selection on r1 and r2
  // leftcol and rightcol point to the column indexes in the BShareTable
  // column 2 is at position 4 and column 3 is at position 6
  // because the BShareTable contains two share-columns per data column
  Predicate_B select_p = {EQ, NULL, NULL, 4, 6};
  select_b(in1, select_p, s_res1);
  select_b(in2, select_p, s_res2);

  // exchange shares of selection results
  exchange_shares_array(s_res1, s2_res1, ROWS1);
  exchange_shares_array(s_res2, s2_res2, ROWS2);

  // 2. Join r1, r2 on id
  Predicate_B join_p = {EQ, NULL, NULL, 0, 0};
  int batch_size = 2;
  BShare result_b[ROWS2*batch_size], remote[ROWS2*batch_size]; // results for one batch
  AShare count = 0;
  AShare converted[ROWS2*batch_size];

  // 3. Logical and of selection results and join result
  // call join_batch for every row of the left input
  for (int i=0, bid=0; i<ROWS1; i+=batch_size, bid++) {
    join_b_batch(&in1, &in2, i, i+batch_size, 0, ROWS2, join_p, remote, result_b);

    // get remote join result
    exchange_shares_array(result_b, remote, ROWS2*batch_size);

    // logical and with s_res2
    for (int j=0; j<ROWS2*batch_size; j++) {
      result_b[j] = and_b(result_b[j], remote[j],
                          s_res2[j%ROWS2], s2_res2[j%ROWS2], get_next_rb());
    }
    
    // get remote and-result
    exchange_shares_array(result_b, remote, ROWS2*batch_size);

    // logical and with s_res1[i]
    int k=i;
    for (int j=0; j<ROWS2*batch_size; j++) {
      result_b[j] = and_b(result_b[j], remote[j],
                          s_res1[k], s2_res1[k], get_next_rb());
      if ((j+1)%ROWS2==0)
        k++;
    }

    // 4. Convert result to arithmetic shares and count selected rows
    convert_single_bit_array(result_b, &ra[bid*ROWS2*batch_size], &rb[bid*ROWS2*batch_size], ROWS2*batch_size, converted);

    // local summing of counts
    for (int k=0; k<ROWS2*batch_size; k++) {
      count += converted[k];
    }
  }
   

  /** ============================================================= *
   *  ASSERT RESULTS
   *  ============================================================= */

  // ASSERT SELECTION
  Data s1_out[ROWS1], s2_out[ROWS2];
  open_b_array(s_res1, ROWS1, s1_out);
  open_b_array(s_res2, ROWS2, s2_out);

  if (rank == 0) {
    for (int i=0; i<ROWS1; i++) {
        if (i!=1 && i!=5) {
          assert(s1_out[i] == 1);
        }
        else {
          assert(s1_out[i] == 0);
        }
    }
    for (int i=0; i<ROWS2; i++) {
        if (i!=1) {
          assert(s2_out[i] == 1);
        }
        else {
          assert(s2_out[i] == 0);
        }
    }
    printf("RELATIONAL (SELECTIONS): OK.\n");
  }
  // ASSERT SELECTION - END

  // ASSERT FINAL RESULT
  Data c_out;
  c_out = open_a(count);
  if (rank == 0) {
    assert(c_out == 2);
    printf("TEST RELATIONAL: OK.\n");
  }

  // ASSERT FINAL RESULT - END

  // tear down communication
  MPI_Finalize();
  return 0;
}
