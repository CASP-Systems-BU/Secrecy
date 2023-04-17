#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS1 3
#define COLS2 1

/**
 * Evaluates the performance of Q1 (comorbidity).
 **/

int main(int argc, char** argv) {

  const long ROWS1 = 8; // input1 size
  const long ROWS2 = 8; // input2 size

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS1, 2*COLS1, 1}; // {pid, diag, cnt}
  allocate_bool_shares_table(&t1);
  BShareTable t2 = {-1, rank, ROWS2, 2*COLS2, 2};
  allocate_bool_shares_table(&t2);

  BShare *rb = (BShare *) malloc(ROWS1*sizeof(BShare));
  AShare *ra = (AShare *) malloc(ROWS1*sizeof(AShare));
  // We need 4 + 6 + 7 (total_comparisons) + 8 (last round equalities) = 25 random numbers 
  BShare *rand_b = (BShare *) malloc(25*sizeof(BShare));
  AShare *rand_a = (AShare *) malloc(25*sizeof(AShare));

  // initialize rand bits for conversion (all equal to 0)
  for (int i=0; i<ROWS1; i++) {
    ra[i] = (unsigned int) 0;
    rb[i] = (unsigned int) 0;
  }

  // initialize rand bits for group-by (all equal to 0)
  for (int i=0; i<25; i++) {
    rand_a[i] = (unsigned int) 0;
    rand_b[i] = (unsigned int) 0;
  }

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data in1[ROWS1][COLS1] = {{0, 10, 0}, {1, 10, 0}, {2, 11, 0}, {3, 9, 0},
                             {4, 7, 0}, {5, 8, 0}, {6, 8, 0}, {7, 8, 0}};
    Data in2[ROWS2][COLS2] = {{0}, {1}, {3}, {5}, {6}, {7}, {30}, {30}};

    Data ** c1 = allocate_2D_data_table(ROWS1, COLS1);
    for (int i=0;i<ROWS1;i++){
      for(int j=0;j<COLS1;j++){
        c1[i][j] = in1[i][j];
      }
    }
    Data ** c2 = allocate_2D_data_table(ROWS2, COLS2);
    for (int i=0;i<ROWS2;i++){
      for(int j=0;j<COLS2;j++){
        c2[i][j] = in2[i][j];
      }
    }
    Table r1 = {-1, ROWS1, COLS1, c1};
    Table r2 = {-1, ROWS2, COLS2, c2};

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

  // start timer
  gettimeofday(&begin, 0);

  // STEP 1: SORT t1 on diag (att=2)
  #if DEBUG
    if (rank==0) {
      printf("Sorting.\n");
    }
  #endif
  unsigned int att_index[1] = {2};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS1/2);

  // STEP 2: IN
  #if DEBUG
    if (rank==0) {
      printf("IN.\n");
    }
  #endif
  BShare *in_res = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(in_res!=NULL);
  in(&t1, &t2, 0, 0, in_res, ROWS1);

  // STEP 3: GROUP-BY-COUNT on diag (att=2)
  #if DEBUG
    if (rank==0) {
      printf("Conversion.\n");
    }
  #endif
  // a. get arithmetic shares of selected bits
  AShare *in_res_a = (AShare *) malloc(ROWS1*sizeof(AShare));
  assert(in_res_a!=NULL);
  convert_single_bit_array(in_res, ra, rb, ROWS1, in_res_a);
  #if DEBUG
    if (rank==0) {
      printf("Group-by.\n");
    }
  #endif
  unsigned key_indices[1] = {2};
  group_by_count_sel_odd_even(&t1, key_indices, 1, 3, in_res, in_res_a, rand_b, rand_a);
  free(rand_a); free(rand_b);

  // reuse ra, rb arrays for exchange of arithmetic counts
  exchange_a_shares_array(in_res_a, ra, ROWS1);

  // STEP 4: sort group's output on count
  // reuse rb, in_res for result of conversion to binary
  #if DEBUG
    if (rank==0) {
      printf("Conversion.\n");
    }
  #endif
  convert_a_to_b_array(in_res_a, ra, rb, in_res, ROWS1);
  free(in_res_a); free(ra);

  // order by count
  // first copy the binary counter to the last column of t1
  for (int i=0; i<ROWS1; i++) {
    t1.content[i][4] = rb[i]; // local share of count
    t1.content[i][5] = in_res[i]; // remote share of count
  }
  free(rb); free(in_res);

  // STEP 5: Sort by cnt
  #if DEBUG
    if (rank==0) {
      printf("Sorting.\n");
    }
  #endif
  att_index[0] = 4;
  asc[0] = 0;
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS1/2);

  // There should be three lines in the result with this order:
  // 8,3    (currently at index 5)
  // 10, 2  (currently at index 6)
  // 9, 1   (currently at index 7)
  // The rest of the lines in the result are garbage (due to multiplexing)
  Data result[ROWS1][2];
  // Open first 8 elements
  for (int i=0; i<ROWS1; i++) {
    result[i][0] = open_b(t1.content[i][2]); // diag
    result[i][1] = open_b(t1.content[i][4]); // count
  }

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    assert(result[5][0]==8);
    assert(result[5][1]==3);
    assert(result[6][0]==10);
    assert(result[6][1]==2);
    assert(result[7][0]==9);
    assert(result[7][1]==1);
    #if DEBUG
      for (int i=0; i<ROWS1; i++) {
        printf("[%d] (diag, cnt) = %lld, %lld\n", i, result[i][0], result[i][1]);
      }
    #endif
  }

  if (rank == 0) {
    printf("TEST Q1: OK.\n");
  }

  free(t1.content); free(t2.content);

  // tear down communication
  MPI_Finalize();
  return 0;
}
