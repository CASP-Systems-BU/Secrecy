#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS 6  // pid, time, time+15, time+56, diag-cdiff, cdiff-diag

/**
 * Evaluates the performance of Q2 (rec. cdiff).
 **/

int main(int argc, char** argv) {

  if (argc < 2) {
    printf("\n\nUsage: %s <NUM_ROWS> \n\n", argv[0]);
    return -1;
  }

  // initialize communication
  init(argc, argv);

  const int ROWS = atol(argv[1]); // input1 size

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  // Diagnosis(pid, time, time+15, time+56, diag-cdif, cdif-diag)
  BShareTable t1 = {-1, rank, ROWS, 2*COLS, 1};
  allocate_bool_shares_table(&t1);

  if (rank == 0) { //P1
    // Initialize input data and shares
    Table r1;
    generate_random_table(&r1, ROWS, COLS);

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS, 2*COLS, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS, 2*COLS, 1};
    allocate_bool_shares_table(&t13);

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r1, &t1, &t12, &t13);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS * 2 * COLS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS * 2 * COLS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1.content);
    free(t12.content);
    free(t13.content);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&(t1.content[0][0]), ROWS * 2 * COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&(t1.content[0][0]), ROWS * 2 * COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  // start timer
  gettimeofday(&begin, 0);

  // STEP 1: Order by pid, time
  #if DEBUG
    if (rank==0) {
      printf("Order-by.\n");
    }
  #endif
  unsigned int att_index[2] = {0,2};
  bool asc[2] = {1,1};
  bitonic_sort_batch(&t1, att_index, 2, asc, ROWS/2);

  // STEP 2: Apply selection predicate 'diag=cdiff'
  #if DEBUG
    if (rank==0) {
      printf("1st selection.\n");
    }
  #endif
  BShare *sel = (BShare *) malloc(ROWS*sizeof(BShare));
  assert(sel!=NULL);
  BShare *rem_sel = (BShare *) malloc(ROWS*sizeof(BShare));
  assert(rem_sel!=NULL);
  // Diagnosis(pid, time, time+15, time+56, diag-cdiff, cdiff-diag)
  Predicate_B p = {EQ, NULL, NULL, 8, 10};
  // Apply 'diag=cdiff'
  select_b(t1, p, sel);

  exchange_shares_array(sel, rem_sel, ROWS);

  // Copy selection bits to 'diag-cdiff' column
  for (int i=0; i<ROWS; i++) {
    t1.content[i][8] = sel[i];
    t1.content[i][9] = rem_sel[i];
  }

  free(sel); free(rem_sel);

  // STEP 3: Apply DISTINCT pid and inequality predicates
  #if DEBUG
    if (rank==0) {
      printf("Distict.\n");
    }
  #endif
  BitShare* d = (BitShare *) malloc(ROWS*sizeof(BitShare));
  assert(d!=NULL);
  BitShare* rem_d = (BitShare *) malloc(ROWS*sizeof(BitShare));
  assert(rem_d!=NULL);
  unsigned key_indices[1] = {0};
  distinct_batch(&t1, key_indices, 1, d, t1.numRows-1);

  exchange_bit_shares_array(d, rem_d, ROWS);

  // Evaluate geq1_i AND geq2_i AND pid1==pid2 AND s_i
  BShare mask=1;
  BitShare *s = (BitShare *) malloc(ROWS*sizeof(BitShare));
  assert(s!=NULL);
  BitShare *rs = (BitShare *) malloc(ROWS*sizeof(BitShare));
  assert(rs!=NULL);
  // Evaluate pid1==pid2 AND s_i AND s_{i+1}
  for (int i=0; i<ROWS-1; i++) {
    // pid[i]==pid[i+1] <=> NOT distinct[i+1]
    s[i] = and_b(d[i+1] ^ mask, rem_d[i+1] ^ mask,
                 t1.content[i][8] ^ mask, t1.content[i][9],
                 get_next_rb())
                & mask;
  }
  s[ROWS-1] = 0;  // Last row is never included in the result
  free(d); free(rem_d);
  
  exchange_bit_shares_array(s, rs, ROWS);

  for (int i=0; i<ROWS-1; i++) {
    // pid[i]==pid[i+1] <=> NOT distinct[i+1]
    s[i] = and_b(s[i], rs[i],
                 t1.content[i + 1][8], t1.content[i + 1][9],
                 get_next_rb())
                & mask;
  }
  exchange_bit_shares_array(s, rs, ROWS);

  // Evaluate second inequality
  BitShare* geq2 = (BitShare *) malloc(t1.numRows*sizeof(BitShare));
  assert(geq2!=NULL);
  BitShare* rem_geq2 = (BitShare *) malloc(t1.numRows*sizeof(BitShare));
  assert(rem_geq2!=NULL);
  adjacent_geq(&t1, 6, 2, geq2, t1.numRows-1, 0);

  exchange_bit_shares_array(geq2, rem_geq2, ROWS);

  // Evaluate geq2_i AND s_i
  for (int i=0; i<ROWS-1; i++) {
    s[i] = and_b(geq2[i], rem_geq2[i], s[i], rs[i], get_next_rb())
                & mask;
  }
  free(geq2); free(rem_geq2);
  exchange_bit_shares_array(s, rs, ROWS);

  // Evaluate first inequality
  BitShare* geq1 = (BitShare *) malloc(ROWS*sizeof(BitShare));
  assert(geq1!=NULL);
  BitShare* rem_geq1 = (BitShare *) malloc(ROWS*sizeof(BitShare));
  assert(rem_geq1!=NULL);
  adjacent_geq(&t1, 2, 4, geq1, t1.numRows-1, 1);

  exchange_bit_shares_array(geq1, rem_geq1, ROWS);
  
  // Evaluate geq1_i AND s_i
  for (int i=0; i<ROWS-1; i++) {
    // Evaluate geq2_i AND s_i
    s[i] = and_b(geq1[i], rem_geq1[i], s[i], rs[i], get_next_rb())
                & mask;
  }
  free(geq1); free(rem_geq1);
  exchange_bit_shares_array(s, rs, ROWS);

  // Update bits in shares table
  for (int i=0; i<ROWS; i++) {
    t1.content[i][8] = s[i];
    t1.content[i][9] = rs[i];
  }

  // STEP 4: Sort diagnosis on s_i, pid
  #if DEBUG
    if (rank==0) {
      printf("2nd Order-by.\n");
    }
  #endif
  att_index[0] = 8;
  att_index[1] = 0;
  bitonic_sort_batch(&t1, att_index, 2, asc, ROWS/2);

   // STEP 5: Scan diagnosis and mask i-th row iff NOT(s_i AND pid[i-1]!=pid[i])
  // Check pid[i]==pid[i+1]
  distinct_batch(&t1, key_indices, 1, s, t1.numRows-1);
  exchange_bit_shares_array(s, rs, ROWS);
  // Evaluate s_i AND NOT(pid[i]==pid[i-1])
  for (int i=0; i<ROWS; i++) {
    s[i] = and_b(t1.content[i][8], t1.content[i][9],
                 s[i] ^ mask, rs[i] ^ mask, get_next_rb())
                & mask;
  }
  exchange_bit_shares_array(s, rs, ROWS);
  BShare b1, b2;
  BShare max=0xFFFFFFFFFFFFFFFF;
  // We only need to multiplex the pid
  BShare *att = (BShare *) malloc(ROWS*sizeof(BShare));
  assert(att!=NULL);
  for (int i=0; i<ROWS; i++) {
    s[i] ^= mask;
    rs[i] ^= mask;
    b1 = - (BShare) s[i];
    b2 = - (BShare) rs[i];
    // Compute pid = b*max + (1-b)*pid
    att[i] = and_b(b1, b2, max, max, get_next_rb());
    att[i] ^= and_b(~b1, ~b2, t1.content[i][0], t1.content[i][1],
                    get_next_rb());
  }

  free(s); free(rs);

  // OPEN diagnosis
  Data *result = (Data *) malloc(ROWS*sizeof(Data));
  assert(result!=NULL);
  open_b_array(att, ROWS, result);

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("\tQ2\t%d\t%.3f\n", ROWS, elapsed);
  }

  free(att);
  free(result); free(t1.content);

  // tear down communication
  MPI_Finalize();
  return 0;
}
