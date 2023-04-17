#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define ROWS1 8
#define COLS1 4

int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // Selected bits
  BShare r[ROWS1] = {1,1,1,0,1,1,1,0};

  // r1(user id1, user id2, score, score)
  BShare r1s1[ROWS1][COLS1], r1s2[ROWS1][COLS1], r1s3[ROWS1][COLS1];
  BShare rs[ROWS1];

  if (rank == 0) { //P1

    // Initialize input data and shares
    Data r1[ROWS1][COLS1] = {{1, 10, 1, 1}, {1, 10, 3, 3}, {1, 10, 6, 6},
                             {4, 11, 2, 12}, {4, 12, 7, 7},
                             {6, 12, 11, 11}, {6, 12, 4, 4}, {6, 12, 2, 2}};
    init_sharing();

    // generate r1 shares
    for (int i=0; i<ROWS1; i++) {
      for (int j=0; j<COLS1; j++) {
        generate_bool_share(r1[i][j], &r1s1[i][j], &r1s2[i][j], &r1s3[i][j]);
      }
    }

    BShare rs2[ROWS1], rs3[ROWS1];
    for (int i=0; i<ROWS1; i++) {
      generate_bool_share(r[i], &rs[i], &rs2[i], &rs3[i]);
    }

    //Send shares to P2
    MPI_Send(&r1s2[0][0], ROWS1*COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r1s3[0][0], ROWS1*COLS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&r1s3[0][0], ROWS1*COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&r1s1[0][0], ROWS1*COLS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    // Send selected bitshares to P2
    MPI_Send(&rs2, ROWS1, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    // Send selected bit shares to P3
    MPI_Send(&rs3, ROWS1, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else { //P2 or_a P3
    MPI_Recv(&r1s1[0][0], ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&r1s2[0][0], ROWS1*COLS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&rs, ROWS1, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  // Copy selected bits because their modified by group_by_min_max_sel()
  BShare rss[ROWS1];
  for (int i=0; i<ROWS1; i++) {
    rss[i] = rs[i];
  }

  // test group_by_join
  #if DEBUG
    if (rank == 0) {
      printf("\nTesting group_by_min_max_sel\n");
    }
  #endif

  // allocate BShare tables
  BShareTable in1 = {-1, rank, ROWS1, 2*COLS1, 1};
  allocate_bool_shares_table(&in1);

  // copy shares into the BShareTables
  for (int i=0; i<ROWS1; i++) {
      in1.content[i][0] = r1s1[i][0];
      in1.content[i][1] = r1s2[i][0];
      in1.content[i][2] = r1s1[i][1];
      in1.content[i][3] = r1s2[i][1];
      in1.content[i][4] = r1s1[i][2];
      in1.content[i][5] = r1s2[i][2];
      in1.content[i][6] = r1s1[i][3];
      in1.content[i][7] = r1s2[i][3];
  }

  unsigned keys[2] = {0,2};
  // in1 (id1, id2, score, score)
  group_by_min_max_sel(&in1, rs, 4, 6, keys, 2);

  // Open the group-by output
  Data out[ROWS1*COLS1], opened[ROWS1*COLS1];
  for (int i=0; i<in1.numRows; i++) {
    for (int j=0; j<in1.numCols; j+=2) {
      out[i*COLS1+j/2] = in1.content[i][j];
    }
  }
  open_b_array(out, ROWS1*COLS1, opened);

  // assert and print result
  if (rank == 0) {
    #if DEBUG
    printf("Output (open): \n");
    #endif
    Data max = 0xFFFFFFFFFFFFFFFF;
    for (int i=0; i<ROWS1*COLS1; i+=COLS1) {
        // r = {1,1,1,0,1,1,1,0};
        // r1 = {{1, 10, 1, 1}, {1, 10, 3, 3}, {1, 10, 6, 6},
        //       {4, 11, 2, 12}, {4, 12, 7, 7},
        //       {6, 12, 11, 11}, {6, 12, 4, 4}, {6, 12, 2, 2}};
        Data gz[ROWS1][COLS1] = {{max, max, max, max}, {max, max, max, max},
                                 {1, 10, 1, 6}, {max, max, max, max},
                                 {4, 12, 7, 7}, {max, max, max, max},
                                 {max, max, max, max}, {6, 12, 4, 11}};
        #if DEBUG
        printf("%lld %lld %lld %lld\n", opened[i], opened[i+1],
                                        opened[i+2], opened[i+3]);
        #endif
        if (i==0 || i==4 || i==12 || i==20 || i==24) {
          assert(opened[i] == max);
        }
        if (i==8) {
          assert(opened[i+2] == gz[2][2]);
          assert(opened[i+3] == gz[2][3]);
        }
        if (i==16) {
          assert(opened[i+2] == gz[4][2]);
          assert(opened[i+3] == gz[4][3]);
        }
        if (i==28) {
          assert(opened[i+2] == gz[7][2]);
          assert(opened[i+3] == gz[7][3]);
        }
    }
    printf("TEST GROUP-BY-MIN-MAX: OK.\n");
  }

  // test group_by_odd_even
  #if DEBUG
    if (rank == 0) {
      printf("\nTesting group_by_min_max_sel_odd_even\n");
    }
  #endif

  // allocate BShare tables
  BShareTable in2 = {-1, rank, ROWS1, 2*COLS1, 1};
  allocate_bool_shares_table(&in2);

  // copy shares into the BShareTables
  for (int i=0; i<ROWS1; i++) {
      in2.content[i][0] = r1s1[i][0];
      in2.content[i][1] = r1s2[i][0];
      in2.content[i][2] = r1s1[i][1];
      in2.content[i][3] = r1s2[i][1];
      in2.content[i][4] = r1s1[i][2];
      in2.content[i][5] = r1s2[i][2];
      in2.content[i][6] = r1s1[i][3];
      in2.content[i][7] = r1s2[i][3];
  }

  // in2 (id1, id2, score, score)
  group_by_min_max_sel_odd_even(&in2, 3, rss, 4, 6, keys, 2);

  // Open the group-by output
  for (int i=0; i<in2.numRows; i++) {
    for (int j=0; j<in2.numCols; j+=2) {
      out[i*COLS1+j/2] = in2.content[i][j];
    }
  }
  open_b_array(out, ROWS1*COLS1, opened);

  // assert and print result
  if (rank == 0) {
    #if DEBUG
    printf("Output (open): \n");
    #endif
    Data max = 0xFFFFFFFFFFFFFFFF;
    for (int i=0; i<ROWS1*COLS1; i+=COLS1) {
        // r = {1,1,1,0,1,1,1,0};
        // r1 = {{1, 10, 1, 1}, {1, 10, 3, 3}, {1, 10, 6, 6},
        //       {4, 11, 2, 12}, {4, 12, 7, 7},
        //       {6, 12, 11, 11}, {6, 12, 4, 4}, {6, 12, 2, 2}};
        Data gz[ROWS1][COLS1] = {{1, 10, 1, 6}, {max, max, max, max},
                                 {max, max, max, max}, {max, max, max, max},
                                 {4, 12, 7, 7}, {6, 12, 4, 11},
                                 {max, max, max, max}, {max, max, max, max}};
        #if DEBUG
        printf("%lld %lld %lld %lld\n", opened[i], opened[i+1],
                                        opened[i+2], opened[i+3]);
        #endif
        if (i==0) {
           assert(opened[0] == gz[0][0]);
           assert(opened[1] == gz[0][1]);
           assert(opened[2] == gz[0][2]);
           assert(opened[3] == gz[0][3]);
        }
        else if (i==4*4) {
          assert(opened[i] == gz[4][0]);
          assert(opened[i+1] == gz[4][1]);
          assert(opened[i+2] == gz[4][2]);
          assert(opened[i+3] == gz[4][3]);
        }
        else if (i==5*4) {
          assert(opened[i] == gz[5][0]);
          assert(opened[i+1] == gz[5][1]);
          assert(opened[i+2] == gz[5][2]);
          assert(opened[i+3] == gz[5][3]);
        }
        else {
          assert(opened[i] == max);
        }
    }
    printf("TEST GROUP-BY-MIN-MAX (ODD-EVEN): OK.\n");
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
