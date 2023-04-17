#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS1 3

/**
 * Evaluates the correctness of Q3 in the Senate paper (password reuse).
 **/

int main(int argc, char** argv) {

  const long ROWS1 = 8; // input size

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input tables per party
  BShareTable t1 = {-1, rank, ROWS1, 2*COLS1, 1}; // {uid, pwd, cnt}
  allocate_bool_shares_table(&t1);

  if (rank == 0) { //P1
    // Initialize input data and shares
    Data in1[ROWS1][COLS1] = {{0, 1, 1}, {0, 1, 1}, {1, 11, 1}, {2, 9, 1},
                             {2, 7, 1}, {5, 8, 1}, {6, 8, 1}, {1, 11, 1}};

    Data ** c1 = allocate_2D_data_table(ROWS1, COLS1);
    for (int i=0;i<ROWS1;i++){
      for(int j=0;j<COLS1;j++){
        c1[i][j] = in1[i][j];
      }
    }

    Table r1 = {-1, ROWS1, COLS1, c1};

    // t1 Bshare tables for P2, P3 (local to P1)
    BShareTable t12 = {-1, 1, ROWS1, 2*COLS1, 1};
    allocate_bool_shares_table(&t12);
    BShareTable t13 = {-1, 2, ROWS1, 2*COLS1, 1};
    allocate_bool_shares_table(&t13);

    init_sharing();

    // Generate boolean shares for r1
    generate_bool_share_tables(&r1, &t1, &t12, &t13);

    //Send shares to P2
    MPI_Send(&(t12.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&(t13.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // free temp tables
    free(r1.content);
    free(t12.content);
    free(t13.content);

  }
  else if (rank == 1) { //P2
    MPI_Recv(&(t1.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&(t1.content[0][0]), ROWS1 * 2 * COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  struct timeval begin, end;
  long seconds, micro;
  double elapsed;

  // start timer
  gettimeofday(&begin, 0);

  // STEP 1: SORT t1 on (uid,pwd) (attributes 1 and 2)
  #if DEBUG
    if (rank==0) {
      printf("Sorting.\n");
    }
  #endif
  unsigned int att_index[2] = {0,2};
  bool asc[2] = {1,1};
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS1/2);

  // STEP 2: GROUP-BY-COUNT on (uid,pwd)
  #if DEBUG
    if (rank==0) {
      printf("Group-by.\n");
    }
  #endif
  unsigned key_indices[2] = {0,2};
  // Create 'selected' bit shares
  BShare sel_b = rank % 2;
  BShare* selected_b = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(selected_b!=NULL);
  for (int i=0; i<ROWS1; i++)
    selected_b[i] = sel_b;
  // Apply group-by and update counters (last two columns) in place
  group_by_sum_rca_sel_odd_even(&t1, 3, selected_b, key_indices, 2);
  free(selected_b);

  // Step 4: Select tuples with count > 1
  #if DEBUG
    if (rank==0) {
      printf("Selection.\n");
    }
  #endif
  BShare sel_s = get_succ() % 2;
  Predicate_B p = {GC, NULL, NULL, 4, -1, sel_b, sel_s};
  BShare* sel_results = (BShare *) malloc(ROWS1*sizeof(BShare));
  select_b(t1, p, sel_results);
  // Step 5: Mask not_b selected tuples
  #if DEBUG
    if (rank==0) {
      printf("Masking.\n");
    }
  #endif
  mask(&t1, sel_results, ROWS1);
  free(sel_results);
  // Step 6: Shuffle before opening -- TODO: This can be done in O(n)
  #if DEBUG
    if (rank==0) {
      printf("Shuffling.\n");
    }
  #endif
  unsigned atts[1] = {0};
  bool asc2[1] = {1};
  bitonic_sort_batch(&t1, atts, 1, asc2, ROWS1/2);

  // There should be two selected lines in the result with this order:
  // 0,1    (currently at index 0)
  // 1, 11  (currently at index 1)
  // The rest of the lines in the result are garbage (due to multiplexing)
  #if DEBUG
    if (rank==0) {
      printf("Open.\n");
    }
  #endif
  // Tuples are opened in batches of size BATCH_SIZE
  int batch_size = 3;
  Data* result = (Data *) malloc(batch_size*2*sizeof(BShare));
  assert(result!=NULL);
  Data* opened = (Data *) malloc(batch_size*2*sizeof(BShare));
  assert(opened!=NULL);
  int start=0, limit=batch_size, step;
  while (start<t1.numRows) {
    for (int i=start, k=0; i<limit; i++, k+=2) {
      result[k] = t1.content[i][0]; // id
      result[k+1] = t1.content[i][2]; // pwd
    }
    open_b_array(result, batch_size*2, opened);
    // Check if result is correct
    if (start==0 && rank==0) {
      assert(opened[0]==0);
      assert(opened[1]==1);
      assert(opened[2]==1);
      assert(opened[3]==11);
      #if DEBUG
        for (int i=0; i<batch_size*2; i+=2) {
          printf("[%d] (id, pwd) = %lld, %lld\n", i, opened[i], opened[i+1]);
        }
      #endif
    }
    // Next step
    start = limit;
    step = limit + batch_size;
    limit = (step <= t1.numRows ? step : t1.numRows);
  }

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  if (rank == 0) {
    printf("TEST PWD_REUSE_BASELINE: OK.\n");
  }

  free(t1.content); free(result); free(opened);

  // tear down communication
  MPI_Finalize();
  return 0;
}
