#include <stdio.h>
#include <assert.h>
#include <sys/time.h>
#include "exp-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS 6  // Three 'original' columns plus two more columns for 'c-a' and
                // 'a-c' shares, and one more column for predicate evaluation

PRIVATE void materialized_join(BShareTable *input1, BShareTable *input2,
                        int leftcol, int rightcol, BShareTable* result);

/**
 * Evaluates the performance of Q3 (aspirin count).
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

  // OFFLINE PHASE: Generate dummy random numbers
  AShare* ra = (AShare *) malloc(ROWS2*ROWS1*sizeof(AShare));
  assert(ra!=NULL);
  BShare* rb = (BShare *) malloc(ROWS2*ROWS1*sizeof(BShare));
  assert(rb!=NULL);
  for (int i=0; i<ROWS2*ROWS1; i++) {
    ra[i] = 0;
    rb[i] = 0;
  }

  // Start timer
  gettimeofday(&begin, 0);

  // STEP 1: Apply selections
  #if DEBUG
    if (rank==0) {
      printf("Applying selections.\n");
    }
  #endif
  BShare *sel1 = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(sel1!=NULL);
  BShare *rem_sel1 = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(rem_sel1!=NULL);
  BShare *sel2 = (BShare *) malloc(ROWS2*sizeof(BShare));
  assert(sel2!=NULL);
  BShare *rem_sel2 = (BShare *) malloc(ROWS2*sizeof(BShare));
  assert(rem_sel2!=NULL);
  // Diagnosis(pid, time, diag, diag-hd, hd-diag, sel)
  Predicate_B p = {EQ, NULL, NULL, 6, 8};
  // Apply 'd.diag=hd'
  select_b(t1, p, sel1);
  // Medication(pid, time, med, med-aspririn, aspririn-med, sel)
  // Apply 'm.med=aspirin'
  select_b(t2, p, sel2);

  exchange_shares_array(sel1, rem_sel1, ROWS1);
  exchange_shares_array(sel2, rem_sel2, ROWS2);

  // Copy shares to extra table columns
  for (int i=0; i<ROWS1; i++) {
    t1.content[i][10] = sel1[i];
    t1.content[i][11] = rem_sel1[i];
  }
  for (int i=0; i<ROWS2; i++) {
    t2.content[i][10] = sel2[i];
    t2.content[i][11] = rem_sel2[i];
  }
  free(sel1); free(sel2); free(rem_sel1); free(rem_sel2);

  // STEP 2: Joins
  #if DEBUG
    if (rank==0) {
      printf("Applying theta-join.\n");
    }
  #endif
  // 'sel2 AND d2 AND sel1 AND d1 AND m.time >= d.time'
  Predicate_B join_pred = {GEQ, NULL, NULL, 2, 2};          // m.time >= d.time
  long res_len = ROWS2 * ROWS1; // Number of elements in the join result
  BShare* result = (BShare *) malloc(res_len*sizeof(BShare)); // join result
  assert(result!=NULL);
  BShare* remote = (BShare *) malloc(res_len*sizeof(BShare)); // remote shares
  assert(remote!=NULL);
  // Used for counting rows
  AShare count=0;
  BShare mask=1;

  // Join t2 with t1 on m.time >= d.time
  #if DEBUG
    if (rank==0) {
      printf("Inequality join.\n");
    }
  #endif
  join_b_batch(&t2, &t1, 0, ROWS2, 0, ROWS1, join_pred, remote,
                result);

  // Get remote join results
  exchange_shares_array(result, remote, res_len);

  #if DEBUG
    if (rank==0) {
      printf("Equality join.\n");
    }
  #endif

  // allocate join result table
  BShareTable res_table = {-1, rank, ROWS1*ROWS2, 2*COLS+2, 1};
  allocate_bool_shares_table(&res_table);

  // Apply join first and materialize result
  materialized_join(&t2, &t1, 0, 0, &res_table);

  // Do the conjuction
  for (int j=0; j<res_len; j++) {
    result[j] = and_b(result[j], remote[j],
                      res_table.content[j][res_table.numCols - 2],
                      res_table.content[j][res_table.numCols - 1], get_next_rb())
                    & mask;
  }

  free(t1.content);

  // Get remote join results
  exchange_shares_array(result, remote, res_len);

  // STEP 3: And with all selection predicates
  // Apply s_j AND d.time <= m.time
  #if DEBUG
    if (rank==0) {
      printf("1st selection.\n");
    }
  #endif
  for (int j=0; j<res_len; j++) {
    result[j] = and_b(result[j], remote[j],
                      res_table.content[j % ROWS1][10],
                      res_table.content[j % ROWS1][11], get_next_rb())
                    & mask;
  }

  // Get remote shares
  exchange_shares_array(result, remote, res_len);

  // Apply s_i AND s_j AND d.time <= m.time
  #if DEBUG
    if (rank==0) {
      printf("2nd selection.\n");
    }
  #endif

  int k=0;
  for (int j=0; j<res_len; j++) {
    result[j] = and_b(result[j], remote[j],
                      t2.content[k][10], t2.content[k][11], get_next_rb())
                    & mask;
    if ((j+1)%ROWS1==0)
      k++;
  }

  // we don't need original tables anymore
  free(t2.content);

  // Get remote shares
  // result, remote contain the selected bits for all predicates
  exchange_shares_array(result, remote, res_len);

  // Copy selected bits into res_table
  for (int j=0; j<res_len; j++) {
    res_table.content[j][res_table.numCols - 2] = result[j];
    res_table.content[j][res_table.numCols - 1] = remote[j];
  }

  free(result); free(remote);

  // sort res_table on selected bit, pid
  #if DEBUG
    if (rank==0) {
      printf("Sorting.\n");
    }
  #endif
  
  // Sort on s_i, pid, time (ASC)
  unsigned int att_ind = res_table.numCols-2;
  unsigned sort_att[2] = {att_ind,0};
  bool asc[2] = {false,true};
  bitonic_sort_batch(&res_table, sort_att, 2, asc, res_len/2);

// Apply DISTINCT pid
  #if DEBUG
    if (rank==0) {
      printf("Applying DISTINCT.\n");
    }
  #endif
  BitShare* d = (BitShare *) malloc(res_len*sizeof(BitShare));
  assert(d!=NULL);
  BShare* rem_d = (BShare *) malloc(res_len*sizeof(BShare));
  assert(rem_d!=NULL);

  unsigned key_indices[1] = {0};
  distinct_batch(&res_table, key_indices, 1, d, res_len-1);

  BShare* dshare = (BShare *) malloc(res_len*sizeof(BShare));
  assert(dshare!=NULL);
  for (int i=0; i<res_len; i++) {
    dshare[i] = (BShare)d[i];
  }

  free(d);
  free(res_table.content);

  // Get remote shares of distinct bits
  exchange_shares_array(dshare, rem_d, res_len);

  // Convert result to arithmetic shares and count selected rows
  #if DEBUG
    if (rank==0) {
      printf("Conversion.\n");
    }
  #endif

  AShare* converted = (AShare *) malloc(res_len*sizeof(AShare));
  assert(converted!=NULL);
  convert_single_bit_array(dshare, ra, rb, res_len, converted);

  // Update local counter
  #if DEBUG
    if (rank==0) {
      printf("Counting.\n");
    }
  #endif
  for (int k=0; k<res_len; k++) {
    count += converted[k];
  }
  free(converted);
  
  // Free memory
  free(rem_d);
  free(ra); free(rb);

  // OPEN COUNT RESULT
  Data o_count = open_a(count);
  #if DEBUG
  if (rank == 0) {
    printf("COUNT: %lld\n", o_count);
  }
  #endif

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("%d\tQ3-BASELINE\t%d\t%d\t%.3f\n",
            COLS, ROWS1, ROWS2, elapsed);
  }

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
    for (int i=0; i<result->numRows; i++) {
      temp_local[i] = result->content[i][res_index];
    }
    exchange_shares_array(temp_local, temp_remote, result->numRows);
      // copy exchanged result back to remote column
    for (int i=0; i<result->numRows; i++) {
      result->content[i][res_index + 1] = temp_remote[i];
    }
  }
  free(temp_local); free(temp_remote);
}