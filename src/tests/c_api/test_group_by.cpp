#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 0
#define SHARE_TAG 193
#define ROWS 10

int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // 'Selected' bits
  Data r[10] = {1, 0, 1, 0, 1, 0, 1, 1, 1, 1};

  AShare converted[10];
  AShare counters[10], remote_counters[10];
  BShare rs[10], rb[10];
  AShare ra[10];

  // We do not need 4*(10-1) = 36 random numbers.
  // We only need 2(10-1) = 18.
  BShare rand_b[36];
  BShare rand_a[36];

  BShare zs1[10][2], zs2[10][2], zs3[10][2];

  // Initialize input data and shares

  Data z[10][2] = {{1, 42}, {1, 42}, {2, 42}, {3, 42}, {15, 42}, {15, 43},
                  {15, 44}, {17, 1}, {18, 1}, {18, 1}};

  if (rank == 0) { //P1

    init_sharing();

    // generate z shares
    for (int i=0; i<10; i++) {
        for (int j=0; j<2; j++) {
            generate_bool_share(z[i][j], &zs1[i][j], &zs2[i][j], &zs3[i][j]);
        }
    }

    BShare rand_b2[36], rand_b3[36];
    AShare rand_a2[36], rand_a3[36];

    BShare rs2[10], rs3[10], rb2[10], rb3[10];
    AShare ra2[10], ra3[10];

    for (int i=0; i<10; i++) {
      generate_bool_share(r[i], &rs[i], &rs2[i], &rs3[i]);
    }

    //Send shares to P2
    MPI_Send(&zs2[0][0], 10*2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs3[0][0], 10*2, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    MPI_Send(&rs2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&zs3[0][0], 10*2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs1[0][0], 10*2, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    MPI_Send(&rs3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // Generate random bits and corresponding shares
    generate_rand_bit_shares(rb, ra, rb2, ra2, rb3, ra3, ROWS);

    generate_rand_bit_shares(rand_b, rand_a, rand_b2,
                             rand_a2, rand_b3, rand_a3, 36);

    // Send random bit shares
    // Send shares to P2
    MPI_Send(&rb2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ra2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rand_b2, 36, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rand_a2, 36, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    // Send shares to P3
    MPI_Send(&rb3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ra3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rand_b3, 36, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rand_a3, 36, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else { //P2 and P3
    MPI_Recv(&zs1[0][0], 10*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&zs2[0][0], 10*2, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&rs, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&rb, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ra, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&rand_b, 36, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&rand_a, 36, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  for (int i=0; i<ROWS; i++) {
    converted[i] = convert_single_bit(rs[i], ra[i], rb[i]);
  }

  // Copy original array of 'selected' bits as it is modified in place
  Data rs4[10];
  memcpy(&rs4, &rs, 10*sizeof(BShare));

  // test group_by_count
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing group_by_count\n", rank);
    }
  #endif

  BShareTable t = {-1, rank, 10, 2*2, 1};
  allocate_bool_shares_table(&t);
  // copy shares into the BShareTables
  for (int i=0; i<10; i++) {
      t.content[i][0] = zs1[i][0];
      t.content[i][1] = zs2[i][0];
      t.content[i][2] = zs1[i][1];
      t.content[i][3] = zs2[i][1];
  }

  // sort in place
  unsigned key_indexes[1] = {0};
  group_by_count(&t, key_indexes, 1, rs, converted, rand_b, rand_a);

  // update share arrays
  BShare zs[20];
  for (int i=0; i<10; i++) {
    zs[2*i] = t.content[i][0];
    zs[2*i+1] = t.content[i][2];
  }

  // reveal the result
  Data out[20];
  open_b_array(zs, 20, out);
  Data c_out[10];
  open_array(converted, 10, c_out);

  if (rank==0) {
    Data max = 0xFFFFFFFFFFFFFFFF;
    // Data r[10] = {1, 0, 1, 0, 1, 0, 1, 1, 1, 1};
    // Data z[10][2] = {{1, 42}, {1, 42}, {2, 42}, {3, 42}, {15, 42}, {15, 43},
    //                 {15, 44}, {17, 1}, {18, 1}, {18, 1}};
    Data g_z[10][2] = {{max, max}, {1, 42}, {2, 42}, {max, max}, {max, max},
                       {max, max}, {15, 44}, {17, 1}, {max, max}, {18, 1}};
    #if DEBUG
      printf("[%d] Grouped array (open):\n", rank);
    #endif
    for (int i=0; i<20; i+=2) {
      #if DEBUG
        printf("(%lld %lld) %lld \t", out[i], out[i+1], c_out[i/2]);
      #endif
      if (i==0 || i==2*3 || i==2*4 || i==2*5 || i==2*8){
        assert(out[i]==max);
      }
      else {
        assert(out[i]==g_z[i/2][0]);
        if (i==2*1 || i==2*2 || i==2*7) {
          assert(c_out[i/2]==1);
        }
        else if (i==2*6 || i==2*9) {
          assert(c_out[i/2]==2);
        }
      }
    }
    printf("TEST GROUP_BY: OK.\n");
  }

  // test group_by_count
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing group_by_count_micro\n", rank);
    }
  #endif

  int succ_rank = get_succ();

  BShareTable t2 = {-1, rank, 10, 2*2, 1};
  allocate_bool_shares_table(&t2);
  // copy shares into the BShareTables
  for (int i=0; i<10; i++) {
      t2.content[i][0] = zs1[i][0];
      t2.content[i][1] = zs2[i][0];
      t2.content[i][2] = zs1[i][1];
      t2.content[i][3] = zs2[i][1];
  }

  for (int i=0; i<ROWS; i++) {
    counters[i] = rank % 2;
    remote_counters[i] = succ_rank % 2;
  }

  // group_by in place
  unsigned keys[1] = {0};
  group_by_count_micro(&t2, keys, 1, counters, remote_counters, rand_b, rand_a);

  // update share arrays
  BShare zs22[20];
  for (int i=0; i<10; i++) {
    zs22[2*i] = t2.content[i][0];
    zs22[2*i+1] = t2.content[i][2];
  }

  // reveal the result
  Data out2[20];
  open_b_array(zs22, 20, out2);
  Data c_out2[10];
  open_array(counters, 10, c_out2);

  if (rank==0) {
    Data max = 0xFFFFFFFFFFFFFFFF;
    // Data z[10][2] = {{1, 42}, {1, 42}, {2, 42}, {3, 42}, {15, 42}, {15, 43},
    //                 {15, 44}, {17, 1}, {18, 1}, {18, 1}};
    Data g_z2[10][2] = {{max, max}, {1, 42}, {2, 42}, {3, 42}, {max, max},
                        {max, max}, {15, 44}, {17, 1}, {max, max}, {18, 1}};
    #if DEBUG
      printf("[%d] Grouped array (open):\n", rank);
    #endif
    for (int i=0; i<20; i+=2) {
      #if DEBUG
        printf("(%lld %lld) %lld \t", out2[i], out2[i+1], c_out2[i/2]);
      #endif
      if (i==0 || i==2*4 || i==2*5 || i==2*8){
        assert(out2[i]==max);
      }
      else {
        assert(out2[i]==g_z2[i/2][0]);
        if (i==2*2 || i==2*3 || i==2*7) {
          assert(c_out2[i/2]==1);
        }
        else if (i==2*1 || i==2*9) {
          assert(c_out2[i/2]==2);
        }
        else if (i==2*6) {
          assert(c_out2[i/2]==3);
        }
      }
    }
    printf("TEST GROUP_BY (MICRO): OK.\n");
  }

  // Test group-by-count using RCA
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing group_by_sum_rca\n", rank);
    }
  #endif

  BShareTable t3 = {-1, rank, 10, 2*2, 1};
  allocate_bool_shares_table(&t3);
  // copy shares into the BShareTables
  for (int i=0; i<10; i++) {
      t3.content[i][0] = zs1[i][0];
      t3.content[i][1] = zs2[i][0];
      t3.content[i][2] = zs1[i][1];
      t3.content[i][3] = zs2[i][1];
  }

  unsigned key_indices[1] = {0};
  // group_by in place
  group_by_sum_rca(&t3, key_indices, 1);

  // update share arrays
  BShare zs32[20];
  for (int i=0; i<10; i++) {
    zs32[2*i] = t3.content[i][0];
    zs32[2*i+1] = t3.content[i][2];
  }

  // reveal the result
  Data out3[20];
  open_b_array(zs32, 20, out3);

  if (rank==0) {
    Data max = 0xFFFFFFFFFFFFFFFF;
    // Data z[10][2] = {{1, 42}, {1, 42}, {2, 42}, {3, 42}, {15, 42}, {15, 43},
    //                 {15, 44}, {17, 1}, {18, 1}, {18, 1}};
    Data g_z2[10][2] = {{max, max}, {1, 42}, {2, 42}, {3, 42}, {max, max},
                        {max, max}, {15, 44}, {17, 1}, {max, max}, {18, 1}};
    #if DEBUG
      printf("[%d] Grouped array (open):\n", rank);
    #endif
    for (int i=0; i<20; i+=2) {
      #if DEBUG
        printf("%lld %lld \t", out3[i], out3[i+1]);
      #endif
      if (i==0 || i==2*4 || i==2*5 || i==2*8){
        assert(out3[i]==max);
      }
      else {
        assert(out3[i]==g_z2[i/2][0]);
        if (i==2*2 || i==2*3) {
          assert(out3[i+1]==42);
        }
        else if (i==2*7) {
          assert(out3[i+1]==1);
        }
        else if (i==2*1) {
          assert(out3[i+1]==84);
        }
        else if (i==2*9) {
          assert(out3[i+1]==2);
        }
        else if (i==2*6) {
          assert(out3[i+1]==129);
        }
      }
    }
    printf("TEST GROUP_BY (RCA): OK.\n");
  }

  // test group_by_sum_rca_sel
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing group_by_sum_rca_sel\n", rank);
    }
  #endif

  BShareTable t4 = {-1, rank, 10, 2*2, 1};
  allocate_bool_shares_table(&t4);
  // copy shares into the BShareTables
  for (int i=0; i<10; i++) {
      t4.content[i][0] = zs1[i][0];
      t4.content[i][1] = zs2[i][0];
      t4.content[i][2] = zs1[i][1];
      t4.content[i][3] = zs2[i][1];
  }

  // sort in place
  unsigned key_indexes4[1] = {0};
  group_by_sum_rca_sel(&t4, rs4, key_indexes4, 1);

  // update share arrays
  BShare zs4[20];
  for (int i=0; i<10; i++) {
    zs4[2*i] = t4.content[i][0];
    zs4[2*i+1] = t4.content[i][2];
  }

  // reveal the result
  Data out4[20];
  open_b_array(zs4, 20, out4);

  if (rank==0) {
    Data max = 0xFFFFFFFFFFFFFFFF;
    // Data r[10] = {1, 0, 1, 0, 1, 0, 1, 1, 1, 1};
    // Data z[10][2] = {{1, 42}, {1, 42}, {2, 42}, {3, 42}, {15, 42}, {15, 43},
    //                 {15, 44}, {17, 1}, {18, 1}, {18, 1}};
    Data g_z[10][2] = {{max, max}, {1, 42}, {2, 42}, {max, max}, {max, max},
                       {max, max}, {15, 86}, {17, 1}, {max, max}, {18, 2}};
    #if DEBUG
      printf("[%d] Grouped array (open):\n", rank);
    #endif
    for (int i=0; i<20; i+=2) {
      #if DEBUG
        printf("%lld %lld \t", out4[i], out4[i+1]);
      #endif
      if (i==0 || i==2*3 || i==2*4 || i==2*5 || i==2*8){
        assert(out4[i]==max);
      }
      else {
        assert(out4[i]==g_z[i/2][0]);
        if (i==2*1 || i==2*2) {
          assert(out4[i+1]==42);
        }
        else if (i==2*6) {
          assert(out4[i+1]==86);
        }
        else if (i==2*7) {
          assert(out4[i+1]==1);
        }
        else if (i==2*9) {
          assert(out4[i+1]==2);
        }
      }
    }
    printf("TEST GROUP_BY (RCA_SEL): OK.\n");
  }

  // tear down communication
  MPI_Finalize();
  return 0;
}
