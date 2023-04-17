#include <stdio.h>
#include <assert.h>
#include <sys/time.h>

#include "exp-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define PRIVATE static
#define COLS1 5

/**
 * Evaluates the correctness of Q4 in the Senate paper (credit score).
 **/

int main(int argc, char** argv) {

  if (argc < 3) {
    printf("\n\nUsage: %s <NUM_ROWS> <BATCH_SIZE>\n\n", argv[0]);
    return -1;
  }

  const int ROWS1 = atol(argv[1]);      // input size
  const int BATCH_SIZE_ = atol(argv[2]); // batch size

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // The input table t(uid, year-2019, 2019-year, score, score+threshold)
  BShareTable t1 = {-1, rank, ROWS1, 2*COLS1, 1};
  allocate_bool_shares_table(&t1);

  if (rank == 0) { //P1
    // Initialize input data and shares
    Table r1;
    generate_random_table(&r1, ROWS1, COLS1);
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

  // STEP 1: SORT t1 on uid (attribute 1)
  #if DEBUG
    if (rank==0) {
      printf("Sorting.\n");
    }
  #endif
  unsigned int att_index[1] = {0};
  bool asc[1] = {1};
  bitonic_sort_batch(&t1, att_index, 1, asc, ROWS1/2);

  // // stop timer
  // gettimeofday(&end, 0);
  // seconds = end.tv_sec - begin.tv_sec;
  // micro = end.tv_usec - begin.tv_usec;
  // elapsed = seconds + micro*1e-6;
  // // Print timing information
  // if (rank == 0) {
  //   printf("%ld\t%.3f\n", ROWS1, elapsed);
  // }

  // STEP 2: Select year=2019
  #if DEBUG
    if (rank==0) {
      printf("Selection.\n");
    }
  #endif
  Predicate_B p = {EQ, NULL, NULL, 2, 4};
  BShare* selected = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(selected!=NULL);
  select_b(t1, p, selected);

  // // stop timer
  // gettimeofday(&end, 0);
  // seconds = end.tv_sec - begin.tv_sec;
  // micro = end.tv_usec - begin.tv_usec;
  // elapsed = seconds + micro*1e-6;
  // // Print timing information
  // if (rank == 0) {
  //   printf("%ld\t%.3f\n", ROWS1, elapsed);
  // }

  // STEP 3: GROUP-BY-MIN-MAX on (uid)
  #if DEBUG
    if (rank==0) {
      printf("Group-by.\n");
    }
  #endif
  unsigned key_indices[1] = {0};
  group_by_min_max_sel_odd_even(&t1, BATCH_SIZE_, selected, 8, 6, key_indices, 1);
  free(selected);

  // // stop timer
  // gettimeofday(&end, 0);
  // seconds = end.tv_sec - begin.tv_sec;
  // micro = end.tv_usec - begin.tv_usec;
  // elapsed = seconds + micro*1e-6;
  // // Print timing information
  // if (rank == 0) {
  //   printf("%ld\t%.3f\n", ROWS1, elapsed);
  // }

  // STEP 4: Select tuples with max-min>threshold
  #if DEBUG
    if (rank==0) {
      printf("Selection.\n");
    }
  #endif
  Predicate_B p2 = {GR, NULL, NULL, 6, 8, -1, -1};
  BShare* sel_results = (BShare *) malloc(ROWS1*sizeof(BShare));
  assert(sel_results!=NULL);
  select_b(t1, p2, sel_results);


  // // stop timer
  // gettimeofday(&end, 0);
  // seconds = end.tv_sec - begin.tv_sec;
  // micro = end.tv_usec - begin.tv_usec;
  // elapsed = seconds + micro*1e-6;
  // // Print timing information
  // if (rank == 0) {
  //   printf("%ld\t%.3f\n", ROWS1, elapsed);
  // }

  // Step 5: Mask not selected tuples
  #if DEBUG
    if (rank==0) {
      printf("Masking.\n");
    }
  #endif
  mask(&t1, sel_results, ROWS1);
  free(sel_results);

  // // stop timer
  // gettimeofday(&end, 0);
  // seconds = end.tv_sec - begin.tv_sec;
  // micro = end.tv_usec - begin.tv_usec;
  // elapsed = seconds + micro*1e-6;
  // // Print timing information
  // if (rank == 0) {
  //   printf("%ld\t%.3f\n", ROWS1, elapsed);
  // }

  // Step 6: Shuffle before opening -- TODO: This can be done in O(n)
  #if DEBUG
    if (rank==0) {
      printf("Shuffling.\n");
    }
  #endif
  unsigned atts[1] = {0};
  bool asc2[1] = {1};
  bitonic_sort_batch(&t1, atts, 1, asc2, ROWS1/2);

  // // stop timer
  // gettimeofday(&end, 0);
  // seconds = end.tv_sec - begin.tv_sec;
  // micro = end.tv_usec - begin.tv_usec;
  // elapsed = seconds + micro*1e-6;
  // // Print timing information
  // if (rank == 0) {
  //   printf("%ld\t%.3f\n", ROWS1, elapsed);
  // }

  // There should be one selected line in the result with this order:
  // 0    (currently at index 0)
  // The rest of the lines in the result are garbage (due to multiplexing)
  #if DEBUG
    if (rank==0) {
      printf("Open.\n");
    }
  #endif
  // Tuples are opened in batches of size BATCH_SIZE_
  Data* result = (Data *) malloc(BATCH_SIZE_*sizeof(BShare));
  assert(result!=NULL);
  Data* opened = (Data *) malloc(BATCH_SIZE_*sizeof(BShare));
  assert(opened!=NULL);
  int start=0, limit=BATCH_SIZE_, step;
  while (start<t1.numRows) {
    for (int i=start, k=0; i<limit; i++, k++) {
      result[k] = t1.content[i][0]; // id
    }
    open_b_array(result, BATCH_SIZE_, opened);
    // Next step
    start = limit;
    step = limit + BATCH_SIZE_;
    limit = (step <= t1.numRows ? step : t1.numRows);
  }

  // stop timer
  gettimeofday(&end, 0);
  seconds = end.tv_sec - begin.tv_sec;
  micro = end.tv_usec - begin.tv_usec;
  elapsed = seconds + micro*1e-6;

  // Print timing information
  if (rank == 0) {
    printf("%d\t%.3f\n", ROWS1, elapsed);
  }

  free(t1.content); free(result); free(opened);

  // tear down communication
  MPI_Finalize();
  return 0;
}
