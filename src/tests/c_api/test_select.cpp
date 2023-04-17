#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define ROWS 10
#define COLS 4
#define EQ_CONST 42

/**
 * SELECT * FROM r WHERE r.1 == 42
**/
int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // r: relation with 10 rows, 4 columns
  BShare rs1[ROWS][COLS], rs2[ROWS][COLS], rs3[ROWS][COLS];

  if (rank == 0) { //P1
    // Initialize input data and shares
    // 3rd column is att-c, 4th is c-att
    Data r[ROWS][COLS] = {{1, 42, 42-EQ_CONST, EQ_CONST-42},
                        {2, 42, 42-EQ_CONST, EQ_CONST-42},
                        {3, 42, 42-EQ_CONST, EQ_CONST-42},
                        {4, 42, 42-EQ_CONST, EQ_CONST-42},
                        {5, 42, 42-EQ_CONST, EQ_CONST-42},
                        {6, 43, 43-EQ_CONST, EQ_CONST-43},
                        {7, 44, 44-EQ_CONST, EQ_CONST-44},
                        {8, 45, 45-EQ_CONST, EQ_CONST-45},
                        {9, 46, 46-EQ_CONST, EQ_CONST-46},
                        {10, 47, 47-EQ_CONST, EQ_CONST-47}
                        };

    init_sharing();

    // generate r shares
    for (int i=0; i<ROWS; i++) {
        for (int j=0; j<COLS; j++) {
            generate_bool_share(r[i][j], &rs1[i][j], &rs2[i][j], &rs3[i][j]);
        }
    }

    //Send shares to P2
    MPI_Send(&rs2[0][0], ROWS*COLS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rs3[0][0], ROWS*COLS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    //Send shares to P3
    MPI_Send(&rs3[0][0], ROWS*COLS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rs1[0][0], ROWS*COLS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else if (rank == 1) { //P2
    MPI_Recv(&rs1[0][0], ROWS*COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&rs2[0][0], ROWS*COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  else { //P3
    MPI_Recv(&rs1[0][0], ROWS*COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&rs2[0][0], ROWS*COLS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  BShare res[ROWS]; // selection result

  // allocate BShare tables
  BShareTable in = {-1, rank, ROWS, 2*COLS, 1};
  allocate_bool_shares_table(&in);

  // copy shares into the BShareTables
  for (int i=0; i<ROWS; i++) {
    for (int j=0; j<COLS; j++) {
        in.content[i][2 * j] = rs1[i][j];
        in.content[i][2 * j + 1] = rs2[i][j];
    }
  }

  // leftcol and rightcol point to the column indexes in the BShareTable
  // column 2 is at position 4 and column 3 is at position 6
  // because the BShareTable contains two share-columns per data column
  Predicate_B p = {EQ, NULL, NULL, 4, 6};

  select_b(in, p, res);

  // reveal the result
  Data out[ROWS]; /** TODO: only rank0 needs to allocate **/
  open_b_array(res, ROWS, out);

  // assert and print result
  if (rank == 0) {
      for (int i=0; i<ROWS; i++) {
          #if DEBUG
            printf("[%d] %lld\t", i, out[i]);
          #endif
          if (i<5) {
              assert(out[i] == 1);
          }
          else {
              assert(out[i] == 0);
          }
      }
      printf("TEST SELECT: OK.\n");
  }

  // leftcol and rightcol point to the column indexes in the BShareTable
  // Create secret-shared constant
  BShare cs1 = (rank==0 ? 43 : (rank==1 ? 0 : 0));
  BShare cs2 = (rank==0 ? 0 : (rank==1 ? 0 : 43));
  Predicate_B p2 = {GC, NULL, NULL, 2, -1, cs1, cs2};

  select_b(in, p2, res);

  // reveal the result
  Data out2[ROWS]; /** TODO: only rank0 needs to allocate **/
  open_b_array(res, ROWS, out2);
  // assert and print result
  if (rank == 0) {
      for (int i=0; i<ROWS; i++) {
          #if DEBUG
            printf("[%d] %lld\t", i, out2[i]);
          #endif
          if (i<=5) {
              assert(out2[i] == 0);
          }
          else {
              assert(out2[i] == 1);
          }
      }
      printf("TEST SELECT_CONST: OK.\n");
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
