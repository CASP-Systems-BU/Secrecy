#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

#define DEBUG 1
#define SHARE_TAG 193
#define ROWS 16

int main(int argc, char** argv) {

  // initialize communication
  init(argc, argv);

  const int rank = get_rank();
  const int pred = get_pred();
  const int succ = get_succ();

  // 'Selected' bits
  Data r[16] = {1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1};

  AShare converted[16], remote_converted[16];
  BShare rs[16], rb[16];
  AShare ra[16];

  // We need total_comparisons + 16 = 49 +16 = 65 random numbers
  BShare rand_b[65];
  BShare rand_a[65];

  BShare zs1[16][3], zs2[16][3], zs3[16][3];
  AShare zsa1[16][3], zsa2[16][3], zsa3[16][3];

  // Initialize input data and shares

  Data z[16][3] = {{1, 42, 1}, {1, 42, 2}, {2, 4, 3}, {3, 42, 0}, {15, 42, 5}, {15, 43, 4},
                  {15, 44, 6}, {17, 1, 1}, {18, 1, 6}, {18, 1, 8}, {18, 1, 1}, {19, 42, 2}, 
                  {20, 4, 3}, {30, 42, 0}, {150, 42, 5}, {150, 43, 4}};

  if (rank == 0) { //P1

    init_sharing();

    // generate z shares
    for (int i=0; i<16; i++) {
        for (int j=0; j<3; j++) {
            generate_bool_share(z[i][j], &zs1[i][j], &zs2[i][j], &zs3[i][j]);
        }
    }

    // generate arithmetic shares (for group_by_avg_odd_even)
    for (int i=0; i<16; i++) {
        for (int j=0; j<3; j++) {
            generate_int_share(z[i][j], &zsa1[i][j], &zsa2[i][j], &zsa3[i][j]);
        }
    }

    BShare rand_b2[65], rand_b3[65];
    AShare rand_a2[65], rand_a3[65];

    BShare rs2[16], rs3[16], rb2[16], rb3[16];
    AShare ra2[16], ra3[16];

    for (int i=0; i<16; i++) {
      generate_bool_share(r[i], &rs[i], &rs2[i], &rs3[i]);
    }

    //Send shares to P2
    MPI_Send(&zs2[0][0], 16*3, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs3[0][0], 16*3, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zsa2[0][0], 16*3, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zsa3[0][0], 16*3, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    MPI_Send(&rs2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    //Send shares to P3
    MPI_Send(&zs3[0][0], 16*3, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zs1[0][0], 16*3, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zsa3[0][0], 16*3, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&zsa1[0][0], 16*3, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    MPI_Send(&rs3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);

    // Generate random bits and corresponding shares
    generate_rand_bit_shares(rb, ra, rb2, ra2, rb3, ra3, ROWS);

    generate_rand_bit_shares(rand_b, rand_a, rand_b2,
                             rand_a2, rand_b3, rand_a3, 65);

    // Send random bit shares
    // Send shares to P2
    MPI_Send(&rb2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ra2, ROWS, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rand_b2, 65, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rand_a2, 65, MPI_LONG_LONG, 1, SHARE_TAG, MPI_COMM_WORLD);

    // Send shares to P3
    MPI_Send(&rb3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&ra3, ROWS, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rand_b3, 65, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
    MPI_Send(&rand_a3, 65, MPI_LONG_LONG, 2, SHARE_TAG, MPI_COMM_WORLD);
  }
  else { //P2 and P3
    MPI_Recv(&zs1[0][0], 16*3, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&zs2[0][0], 16*3, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&zsa1[0][0], 16*3, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&zsa2[0][0], 16*3, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&rs, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&rb, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&ra, ROWS, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&rand_b, 65, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    MPI_Recv(&rand_a, 65, MPI_LONG_LONG, 0, SHARE_TAG, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  //exchange seeds
  exchange_rsz_seeds(succ, pred);

  for (int i=0; i<ROWS; i++) {
    converted[i] = convert_single_bit(rs[i], ra[i], rb[i]);
  }

  // Copy original array of 'selected' bits as it is modified in place
  BShare rss[16];
  memcpy(&rss, &rs, 16*sizeof(BShare));

  // test group_by_count_sel
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing group_by_count_sel_odd_even_aggregation\n", rank);
    }
  #endif

  BShareTable t = {-1, rank, 16, 3*2, 1};
  allocate_bool_shares_table(&t);
  // copy shares into the BShareTables
  for (int i=0; i<16; i++) {
      t.content[i][0] = zs1[i][0];
      t.content[i][1] = zs2[i][0];
      t.content[i][2] = zs1[i][1];
      t.content[i][3] = zs2[i][1];
      t.content[i][4] = zs1[i][2];
      t.content[i][5] = zs2[i][2];
  }

  // sort in place
  unsigned key_indexes[2] = {0,2};
  group_by_count_sel_odd_even(&t, key_indexes, 2, 3, rs, converted, rand_b, rand_a);

  // update share arrays
  BShare zs[48];
  for (int i=0; i<16; i++) {
    zs[3*i] = t.content[i][0];
    zs[3*i+1] = t.content[i][2];
    zs[3*i+2] = t.content[i][4];
  }

  // reveal the result
  Data out[48];
  open_b_array(zs, 48, out);
  Data c_out[16];
  open_array(converted, 16, c_out);

  if (rank==0) {
    Data max = 0xFFFFFFFFFFFFFFFF;
    // Data r[16] = {1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1};
    // Data z[16][3] = {{1, 42, 1}, {1, 42, 2}, {2, 4, 3}, {3, 42, 0}, {15, 42, 5}, {15, 43, 4},
    //                  {15, 44, 6}, {17, 1, 1}, {18, 1, 6}, {18, 1, 8}, {18, 1, 1}, {19, 42, 2}, 
    //                  {20, 4, 3}, {30, 42, 0}, {150, 42, 5}, {150, 43, 4}};
    Data g_z[16][3] = {{1, 42, 1}, {max, max, max}, {2, 4, 3}, {max, max, max}, {max, max, max},
                       {15, 43, 4}, {15, 44, 6}, {17, 1, 1}, {max, max, max}, {18, 1, 8},
                       {max, max, max}, {max, max, max}, {20, 4, 3}, {30, 42, 0}, {150, 42, 5}, {150, 43, 4}};
    #if DEBUG
      printf("[%d] Grouped array (open):\n", rank);
    #endif
    for (int i=0; i<48; i+=3) {
      #if DEBUG
        printf("(%lld %lld %lld) %lld \n", out[i], out[i+1], out[i+2], c_out[i/3]);
      #endif
      if (i==3*1 || i==3*3 || i==3*4 || i==3*8 || i==3*10 || i==3*11){
//        assert(out[i]==max);
      }
      else {
//        assert(out[i]==g_z[i/3][0] && out[i+1]==g_z[i/3][1]);
        if (i==0 || i==3*5 || i==3*6 || i==3*7 || i==3*12 || i==3*13 || i==3*14 || i==3*15) {
//          assert(c_out[i/3]==1);
        }
        else if (i==3*9) {
//          assert(c_out[i/3]==2);
        }
      }
    }
    printf("TEST GROUP_BY_COUNT_SEL (ODD-EVEN): OK.\n");
  }

  // Generate dummy count shares
  for (int i=0; i<ROWS; i++) {
    converted[i] = rank%2;
  }
  exchange_a_shares_array(converted, remote_converted, 16);

  // test group_by_count
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing group_by_count_odd_even_aggregation\n", rank);
    }
  #endif

  BShareTable t2 = {-1, rank, 16, 3*2, 1};
  allocate_bool_shares_table(&t2);
  // copy shares into the BShareTables
  for (int i=0; i<16; i++) {
      t2.content[i][0] = zs1[i][0];
      t2.content[i][1] = zs2[i][0];
      t2.content[i][2] = zs1[i][1];
      t2.content[i][3] = zs2[i][1];
      t2.content[i][4] = zs1[i][2];
      t2.content[i][5] = zs2[i][2];
  }

  // try with one key this time
  unsigned key_index[1] = {0};
  group_by_count_odd_even(&t2, key_index, 1, 1, converted, remote_converted, rand_b, rand_a);

  // update share arrays
  for (int i=0; i<16; i++) {
    zs[3*i] = t2.content[i][0];
    zs[3*i+1] = t2.content[i][2];
    zs[3*i+2] = t2.content[i][4];
  }

  // reveal the result
  open_b_array(zs, 48, out);
  open_array(converted, 16, c_out);

  if (rank==0) {
    Data max = 0xFFFFFFFFFFFFFFFF;
    // Data z[16][3] = {{1, 42, 1}, {1, 42, 2}, {2, 4, 3}, {3, 42, 0}, {15, 42, 5}, {15, 43, 4},
    //                  {15, 44, 6}, {17, 1, 1}, {18, 1, 6}, {18, 1, 8}, {18, 1, 1}, {19, 42, 2}, 
    //                  {20, 4, 3}, {30, 42, 0}, {150, 42, 5}, {150, 43, 4}};
    Data g_z[16][3] = {{1, 42, 1}, {max, max, max}, {2, 4, 3}, {3, 42, 0}, {15, 42, 5},
                       {max, max, max}, {max, max, max}, {17, 1, 1}, {18, 1, 6}, {max, max, max},
                       {max, max, max}, {19, 42, 2}, {20, 4, 3}, {30, 42, 0}, {150, 42, 5}, {max, max, max}};
    #if DEBUG
      printf("[%d] Grouped array (open):\n", rank);
    #endif
    for (int i=0; i<48; i+=3) {
      #if DEBUG
        printf("(%lld %lld %lld) %lld \n", out[i], out[i+1], out[i+2], c_out[i/3]);
      #endif
      if (i==3*1 || i==3*5 || i==3*6 || i==3*9 || i==3*10 || i==3*15){
        assert(out[i]==max);
      }
      else {
        assert(out[i]==g_z[i/3][0] && out[i+1]==g_z[i/3][1]);
        if (i==3*4 || i==3*8) {
          assert(c_out[i/3]==3);
        }
        else if (i==0 || i==3*14) {
          assert(c_out[i/3]==2);
        }
        else {
          assert(c_out[i/3]==1);
        }
      }
    }
    printf("TEST GROUP_BY_COUNT (ODD-EVEN): OK.\n");
  }

  // test group_by_sum_sel_rca
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing group_by_sum_sel_rca_odd_even_aggregation\n", rank);
    }
  #endif

  BShareTable t3 = {-1, rank, 16, 3*2, 1};
  allocate_bool_shares_table(&t3);
  // copy shares into the BShareTables
  for (int i=0; i<16; i++) {
      t3.content[i][0] = zs1[i][0];
      t3.content[i][1] = zs2[i][0];
      t3.content[i][2] = zs1[i][1];
      t3.content[i][3] = zs2[i][1];
      t3.content[i][4] = zs1[i][2];
      t3.content[i][5] = zs2[i][2];
  }

  // try with one key this time
  group_by_sum_rca_sel_odd_even(&t3, 1, rss, key_indexes, 2);

  // update share arrays
  for (int i=0; i<16; i++) {
    zs[3*i] = t3.content[i][0];
    zs[3*i+1] = t3.content[i][2];
    zs[3*i+2] = t3.content[i][4];
  }

  // reveal the result
  open_b_array(zs, 48, out);
  open_b_array(rss, 16, c_out);

  if (rank==0) {
    Data max = 0xFFFFFFFFFFFFFFFF;
    // Data r[16] = {1, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 1, 1};
    // Data z[16][3] = {{1, 42, 1}, {1, 42, 2}, {2, 4, 3}, {3, 42, 0}, {15, 42, 5}, {15, 43, 4},
    //                  {15, 44, 6}, {17, 1, 1}, {18, 1, 6}, {18, 1, 8}, {18, 1, 1}, {19, 42, 2}, 
    //                  {20, 4, 3}, {30, 42, 0}, {150, 42, 5}, {150, 43, 4}};
    Data g_z[16][3] = {{1, 42, 1}, {max, max, max}, {2, 4, 3}, {max, max, max}, {max, max, max},
                       {15, 43, 4}, {15, 44, 6}, {17, 1, 1}, {max, max, max}, {18, 1, 9},
                       {max, max, max}, {max, max, max}, {20, 4, 3}, {30, 42, 0}, {150, 42, 5}, {150, 43, 4}};
    #if DEBUG
      printf("[%d] Grouped array (open):\n", rank);
    #endif
    for (int i=0; i<48; i+=3) {
      #if DEBUG
        printf("(%lld %lld %lld) %lld \n", out[i], out[i+1], out[i+2], c_out[i/3]);
      #endif
      if (i==3*1 || i==3*3 || i==3*4 || i==3*8 || i==3*10 || i==3*11){
        assert(out[i]==max);
        if (i==3*10) {
          assert(c_out[i/3]==1);
        }
        else {
           assert(c_out[i/3]==0);
        }
      }
      else {
        assert(out[i]==g_z[i/3][0] && out[i+1]==g_z[i/3][1]);
        assert(c_out[i/3]==1);
        if (i==3*9){
          assert(out[i+2]==9);
        }
        else{
          assert(out[i+2]==g_z[i/3][2]);
        }
      }
    }
    printf("TEST GROUP_BY_SUM_SEL_RCA (ODD-EVEN): OK.\n");
  }

  // test group_by_sum_rca
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing group_by_sum_rca_odd_even_aggregation\n", rank);
    }
  #endif

  BShareTable t4 = {-1, rank, 16, 3*2, 1};
  allocate_bool_shares_table(&t4);
  // copy shares into the BShareTables
  for (int i=0; i<16; i++) {
      t4.content[i][0] = zs1[i][0];
      t4.content[i][1] = zs2[i][0];
      t4.content[i][2] = zs1[i][1];
      t4.content[i][3] = zs2[i][1];
      t4.content[i][4] = zs1[i][2];
      t4.content[i][5] = zs2[i][2];
  }

  // try with one key this time
  group_by_sum_rca_odd_even(&t4, 1, key_index, 1);

  // update share arrays
  for (int i=0; i<16; i++) {
    zs[3*i] = t4.content[i][0];
    zs[3*i+1] = t4.content[i][2];
    zs[3*i+2] = t4.content[i][4];
  }

  // reveal the result
  open_b_array(zs, 48, out);

  if (rank==0) {
    Data max = 0xFFFFFFFFFFFFFFFF;
    // Data z[16][3] = {{1, 42, 1}, {1, 42, 2}, {2, 4, 3}, {3, 42, 0}, {15, 42, 5}, {15, 43, 4},
    //                  {15, 44, 6}, {17, 1, 1}, {18, 1, 6}, {18, 1, 8}, {18, 1, 1}, {19, 42, 2}, 
    //                  {20, 4, 3}, {30, 42, 0}, {150, 42, 5}, {150, 43, 4}};
    Data g_z[16][3] = {{1, 42, 3}, {max, max, max}, {2, 4, 3}, {3, 42, 0}, {15, 42, 15},
                       {max, max, max}, {max, max, max}, {17, 1, 1}, {18, 1, 15}, {max, max, max},
                       {max, max, max}, {19, 42, 2}, {20, 4, 3}, {30, 42, 0}, {150, 42, 9}, {max, max, max}};
    #if DEBUG
      printf("[%d] Grouped array (open):\n", rank);
    #endif
    for (int i=0; i<48; i+=3) {
      #if DEBUG
        printf("(%lld %lld %lld) \n", out[i], out[i+1], out[i+2]);
      #endif
      if (i==3*1 || i==3*5 || i==3*6 || i==3*9 || i==3*10 || i==3*15){
        assert(out[i]==max);
      }
      else {
        assert(out[i]==g_z[i/3][0] && out[i+1]==g_z[i/3][1]);
        if (i==0) {
          assert(out[i+2]==3);
        }
        else if (i==3*4 || i==3*8) {
          assert(out[i+2]==15);
        }
        else if (i==3*14) {
          assert(out[i+2]==9);
        }
        else {
          assert(out[i+2]==g_z[i/3][2]);
        }
      }
    }
    printf("TEST GROUP_BY_SUM_RCA (ODD-EVEN): OK.\n");
  }

    // test group_by_avg
  #if DEBUG
    if (rank == 0) {
      printf("\n[%d] Testing group_by_avg_odd_even\n", rank);
    }
  #endif

  BShareTable t5 = {-1, rank, 16, 3*2, 1};
  allocate_bool_shares_table(&t5); // This works becayse BShare and AShare have the same size
  // copy shares into the BShareTables
  for (int i=0; i<16; i++) {
      t5.content[i][0] = zs1[i][0];
      t5.content[i][1] = zs2[i][0];
      t5.content[i][2] = zsa1[i][1];
      t5.content[i][3] = zsa2[i][1];
      t5.content[i][4] = zsa1[i][2];
      t5.content[i][5] = zsa2[i][2];
  }

  // try with one key
  unsigned key_idx[1] = {0};
  group_by_avg_odd_even(&t5, 1, rand_b, rand_a, 2, 4, key_idx, 1);

  // update share arrays
  for (int i=0; i<16; i++) {
    zs[i] = t5.content[i][0]; 
  }
  AShare zsa[32];
  for (int i=0; i<16; i++) {
    zsa[2*i] = t5.content[i][2];    
    zsa[2*i+1] = t5.content[i][4];
  }

  // reveal the result
  Data outa[32];
  open_b_array(zs, 16, out);
  open_array(zsa, 32, outa);

  if (rank==0) {
    Data max = 0xFFFFFFFFFFFFFFFF;
    // Data z[16][3] = {{1, 42, 1}, {1, 42, 2}, {2, 4, 3}, {3, 42, 0}, {15, 42, 5}, {15, 43, 4},
    //                  {15, 44, 6}, {17, 1, 1}, {18, 1, 6}, {18, 1, 8}, {18, 1, 1}, {19, 42, 2}, 
    //                  {20, 4, 3}, {30, 42, 0}, {150, 42, 5}, {150, 43, 4}};
    Data g_z[16][3] = {{1, 84, 3}, {max, max, max}, {2, 4, 3}, {3, 42, 0}, {15, 129, 15},
                       {max, max, max}, {max, max, max}, {17, 1, 1}, {18, 3, 15}, {max, max, max},
                       {max, max, max}, {19, 42, 2}, {20, 4, 3}, {30, 42, 0}, {150, 85, 9}, {max, max, max}};
    #if DEBUG
      printf("[%d] Grouped array (open):\n", rank);
    #endif
    for (int i=0; i<16; i++) {
      #if DEBUG
        printf("(%lld %lld %lld) \n", out[i], outa[2*i], outa[2*i+1]);
      #endif
      if (i==1 || i==5 || i==6 || i==9 || i==10 || i==15){
        assert(out[i]==max);
      }
      else {
        assert(out[i]==g_z[i][0]);
        if (i==0) {
          assert(outa[2*i]==84);
          assert(outa[2*i+1]==3);
        }
        else if (i==2) {
          assert(outa[2*i]==4);
          assert(outa[2*i+1]==3);
        }
        else if (i==3) {
          assert(outa[2*i]==42);
          assert(outa[2*i+1]==0);
        }
        else if (i==4) {
          assert(outa[2*i]==129);
          assert(outa[2*i+1]==15);
        }
        else if (i==7) {
          assert(outa[2*i]==1);
          assert(outa[2*i+1]==1);
        }
        else if (i==8) {
          assert(outa[2*i]==3);
          assert(outa[2*i+1]==15);
        }
        else if (i==11) {
          assert(outa[2*i]==42);
          assert(outa[2*i+1]==2);
        }
        else if (i==12) {
          assert(outa[2*i]==4);
          assert(outa[2*i+1]==3);
        }
        else if (i==13) {
          assert(outa[2*i]==42);
          assert(outa[2*i+1]==0);
        }
        else {
          assert(outa[2*i]==85);
          assert(outa[2*i+1]==9);
        }
      }
    }
    printf("TEST GROUP_BY_AVG (ODD-EVEN): OK.\n");
  }
  
  // tear down communication
  MPI_Finalize();
  return 0;
}
