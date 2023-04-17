#include "../../include/core/baseline.h"

#define PRIVATE static

PRIVATE unsigned long long geq_round_a(BShare, BShare, BShare, BShare, int);
PRIVATE unsigned long long gr_round_a(BShare, BShare, BShare, BShare, int);
PRIVATE unsigned long long gr_round_b(BShare, BShare, BShare, BShare, int,
                                      char local[], char remote[]);
PRIVATE unsigned long long gr_round_c(int, int, int, char local[],
                                      char remote[], int levels[],
                                      int *bit_count);
PRIVATE unsigned long long gr_round_c_char(int, int, int, char local[],
                                           char remote[], char levels[],
                                           int *bit_count);
PRIVATE void eq_bulk(int, int, BShare*, BShare*, int);


// array-based boolean equality between a vector and a constant
void eq_b_array_const(BShare *x1, BShare *x2, BShare y1, BShare y2, long len,
                      BShare *res) {

  BShare *res2 = (BShare*) malloc(len*sizeof(BShare)); // remote shares
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);

  // compute bitwise x^y^1
  for (long i=0; i<len; i++) {
    res[i] = x1[i] ^ y1 ^ (~(BShare)0); // local share;
    res2[i] = x2[i] ^ y2 ^ (~(BShare)0); // remote share
  }

  // The result is stored in the (numbits/2) rightmost bits of res, res2 elements
  for (int l=0; l<numlevels; l++) {

    for (long i=0; i<len; i++) {
      res[i] = eq_b_level2(numbits >> l, res[i], res2[i]);
    }

    // Exchange results of logical and, except for the final round
    if (l != numlevels-1) {
      exchange_shares_array(res, res2, len);
    }

  }
  // The local share of the final result is stored in res.
  // We need to call exchange again before using it in a subsequent operation.
  free(res2);
}

// Greater or equal with a secret-shared constant
void geq_batch_const(const BShare *x1, const BShare *x2, BShare y1, BShare y2,
                     int numElements, BitShare *res) {

  int share_length=sizeof(BShare)*8; // The length of BShare in number of bits
  int len = (share_length-1)*sizeof(char) + sizeof(BShare);
  // Local and remote bits per level
  unsigned long long *local_bits = (unsigned long long *) malloc(numElements * sizeof(BShare));
  assert(local_bits!=NULL);
  unsigned long long *remote_bits = (unsigned long long *) malloc(numElements * sizeof(BShare));
  assert(remote_bits!=NULL);

  // For each element, we reserve 1 byte for 63 levels + 1 BShare for last level
  char **local = allocate_2D_byte_array(numElements, len);
  char **remote = allocate_2D_byte_array(numElements, len);

  /** FIRST ROUND **/
  for (int i=0; i<numElements; i++) {
    local_bits[i] = geq_round_a(x1[i], x2[i], y1, y2, share_length);
  }

  // Get the second share of each bit as computed by the other party
  exchange_shares_array_u(local_bits, remote_bits, numElements);

  // Unpack bits and update levels
  for (int i=0; i<numElements; i++) {
    for (int j=0; j<share_length-1; j++) {
      local[i][j] = get_bit_u(local_bits[i], share_length-j-1);
      remote[i][j] = get_bit_u(remote_bits[i], share_length-j-1);
    }
    // Update last-level bits
    unsigned long long l_tmp = get_bit_u(local_bits[i], 0);
    unsigned long long r_tmp = get_bit_u(remote_bits[i], 0);
    memcpy(&local[i][share_length-1], &l_tmp, sizeof(unsigned long long));
    memcpy(&remote[i][share_length-1], &r_tmp, sizeof(unsigned long long));
  }

  /** SECOND ROUND **/
  for (int i=0; i<numElements; i++) {
    local_bits[i] = gr_round_b(x1[i], x2[i], y1, y2, share_length,
                    &local[i][0], &remote[i][0]);
  }

  // Get the second share of each bit as computed by the other party
  exchange_shares_array_u(local_bits, remote_bits, numElements);

  for (int i=0; i<numElements; i++) {
    // Unpack the length/2 MSBs and store them at the last level
    unsigned long long tmp = ( remote_bits[i] >> (share_length/2) );
    memcpy(&remote[i][share_length-1], &tmp, sizeof(unsigned long long));
    // Unpack the rest and update odd levels
    for (int j=1; j<share_length-1; j+=2) {
      remote[i][j] = get_bit_u(remote_bits[i], j/2);
    }
  }

  /** REMAINING ROUNDS **/
  int rounds = (int) log2(share_length/2);
  int **levels = allocate_int_2D_table(numElements, share_length);  // max 'length' levels per pair

  // Initialize level cache
  for (int i=0; i<numElements; i++) {
    for (int j=0; j<share_length; j++) {
      levels[i][j] = -1;
    }
  }

  int bits_left=share_length, bit_count;
  for (int r=1; r<=rounds; r++) {
    bits_left /= 2;
    for (int i=0; i<numElements; i++) {
      bit_count = 0;
      local_bits[i] = gr_round_c(r, bits_left, share_length,
                                 &local[i][0], &remote[i][0],
                                 &levels[i][0], &bit_count);
    }

    // Exchange all bits of the current round and unpack accordingly
    exchange_shares_array_u(local_bits, remote_bits, numElements);

    // Unpack bits of last level
    for (int i = 0; i < numElements; i++) {
      unsigned long long tmp = ( remote_bits[i] >> bit_count );
      memcpy(&remote[i][share_length-1], &tmp, sizeof(BShare));
      // Unpack the rest and reset level cache for next round
      int l=0;
      while ((levels[i][l] >= 0) & (l < share_length)) {
        remote[i][levels[i][l]] = get_bit_u(remote_bits[i], l);
        levels[i][l++] = -1;  // Reset for next round
      }
    }
  }

  free(local_bits); free(remote_bits); free(levels);

  // One bitshare for each greater() comparison
  BitShare mask = 1;
  for (int i=0; i<numElements; i++) {
    res[i] = 0;
    // Do a final XOR of all levels
    for (int j=0; j<share_length-1; j++) {
      res[i] ^= local[i][j];
    }
    // XOR with last level
    res[i] ^= *((unsigned long long*) &local[i][share_length-1]);
    res[i] &= mask;
  }

  free(local); free(remote);
}

// Greater than a secret-shared constant
void greater_batch_const(const BShare *x1, const BShare *x2,
                         const BShare y1, const BShare y2,
                         int numElements, BitShare *res) {

  int share_length=sizeof(BShare)*8; // The length of BShare in number of bits
  int len = (share_length-1)*sizeof(char) + sizeof(BShare);
  // Local and remote bits per level
  unsigned long long *local_bits = (unsigned long long *) malloc(numElements * sizeof(BShare));
  assert(local_bits!=NULL);
  unsigned long long *remote_bits = (unsigned long long *) malloc(numElements * sizeof(BShare));
  assert(remote_bits!=NULL);

  // For each element, we reserve 1 byte for 63 levels + 1 BShare for last level
  char **local = allocate_2D_byte_array(numElements, len);
  char **remote = allocate_2D_byte_array(numElements, len);

  /** FIRST ROUND **/
  for (int i=0; i<numElements; i++) {
    local_bits[i] = gr_round_a(x1[i], x2[i], y1, y2, share_length);
  }

  // Get the second share of each bit as computed by the other party
  exchange_shares_array_u(local_bits, remote_bits, numElements);

  // Unpack bits and update levels
  for (int i=0; i<numElements; i++) {
    for (int j=0; j<share_length-1; j++) {
      local[i][j] = get_bit_u(local_bits[i], share_length-j-1);
      remote[i][j] = get_bit_u(remote_bits[i], share_length-j-1);
    }
    // Update last-level bits
    unsigned long long l_tmp = get_bit_u(local_bits[i], 0);
    unsigned long long r_tmp = get_bit_u(remote_bits[i], 0);
    memcpy(&local[i][share_length-1], &l_tmp, sizeof(unsigned long long));
    memcpy(&remote[i][share_length-1], &r_tmp, sizeof(unsigned long long));
  }

  /** SECOND ROUND **/
  for (int i=0; i<numElements; i++) {
    local_bits[i] = gr_round_b(x1[i], x2[i], y1, y2, share_length,
                    &local[i][0], &remote[i][0]);
  }

  // Get the second share of each bit as computed by the other party
  exchange_shares_array_u(local_bits, remote_bits, numElements);

  for (int i=0; i<numElements; i++) {
    // Unpack the length/2 MSBs and store them at the last level
    unsigned long long tmp = ( remote_bits[i] >> (share_length/2) );
    memcpy(&remote[i][share_length-1], &tmp, sizeof(unsigned long long));
    // Unpack the rest and update odd levels
    for (int j=1; j<share_length-1; j+=2) {
      remote[i][j] = get_bit_u(remote_bits[i], j/2);
    }
  }

  /** REMAINING ROUNDS **/
  int rounds = (int) log2(share_length/2);
  char **levels = allocate_2D_byte_array(numElements, share_length/2);  // max 'length' levels per pair

  // Initialize level cache
  for (int i=0; i<numElements; i++) {
    for (int j=0; j<share_length/2; j++) {
      levels[i][j] = -1;
    }
  }

  int bits_left=share_length, bit_count;
  for (int r=1; r<=rounds; r++) {
    bits_left /= 2;
    for (int i=0; i<numElements; i++) {
      bit_count = 0;
      local_bits[i] = gr_round_c_char(r, bits_left, share_length,
                                 &local[i][0], &remote[i][0],
                                 &levels[i][0], &bit_count);
    }

    // Exchange all bits of the current round and unpack accordingly
    exchange_shares_array_u(local_bits, remote_bits, numElements);

    // Unpack bits of last level
    for (int i = 0; i < numElements; i++) {
      unsigned long long tmp = ( remote_bits[i] >> bit_count );
      memcpy(&remote[i][share_length-1], &tmp, sizeof(BShare));
      // Unpack the rest and reset level cache for next round
      int l=0;
      while ((levels[i][l] >= 0) & (l < share_length/2)) {
        remote[i][(int)levels[i][l]] = get_bit_u(remote_bits[i], l);
        levels[i][l++] = -1;  // Reset for next round
      }
    }
  }

  free(local_bits); free(remote_bits); free(levels);

  // One bitshare for each greater() comparison
  BitShare mask = 1;
  for (int i=0; i<numElements; i++) {
    res[i] = 0;
    // Do a final XOR of all levels
    for (int j=0; j<share_length-1; j++) {
      res[i] ^= local[i][j];
    }
    // XOR with last level
    res[i] ^= *((unsigned long long*) &local[i][share_length-1]);
    res[i] &= mask;
  }

  free(local); free(remote);
}

// Checks if each element in the specified column is greater than the given
// secret-shared constant
void select_greater_batch_const(BShareTable input, int leftcol,
                                BShare cshare1, BShare cshare2,
                                BShare result[]) {

  int share_length=sizeof(BShare)*8; // The length of BShare in number of bits
  int len = (share_length-1)*sizeof(char) + sizeof(BShare);
  int numElements = input.numRows;
  BShare** x = input.content;
  // Local and remote bits per level
  unsigned long long *local_bits = (unsigned long long *) malloc(numElements * sizeof(BShare));
  assert(local_bits!=NULL);
  unsigned long long *remote_bits = (unsigned long long *) malloc(numElements * sizeof(BShare));
  assert(remote_bits!=NULL);

  // For each element, we reserve 1 byte for 63 levels + 1 BShare for last level
  char **local = allocate_2D_byte_array(numElements, len);
  char **remote = allocate_2D_byte_array(numElements, len);

  /** FIRST ROUND **/
  for (int i=0; i<numElements; i++) {
    local_bits[i] = gr_round_a(x[i][leftcol], x[i][leftcol+1], cshare1, cshare2,
                               share_length);
  }

  // Get the second share of each bit as computed by the other party
  exchange_shares_array_u(local_bits, remote_bits, numElements);

  // Unpack bits and update levels
  for (int i=0; i<numElements; i++) {
    for (int j=0; j<share_length-1; j++) {
      local[i][j] = get_bit_u(local_bits[i], share_length-j-1);
      remote[i][j] = get_bit_u(remote_bits[i], share_length-j-1);
    }
    // Update last-level bits
    unsigned long long l_tmp = get_bit_u(local_bits[i], 0);
    unsigned long long r_tmp = get_bit_u(remote_bits[i], 0);
    memcpy(&local[i][share_length-1], &l_tmp, sizeof(unsigned long long));
    memcpy(&remote[i][share_length-1], &r_tmp, sizeof(unsigned long long));
  }

  /** SECOND ROUND **/
  for (int i=0; i<numElements; i++) {
    local_bits[i] = gr_round_b(x[i][leftcol], x[i][leftcol+1], cshare1, cshare2,
                               share_length, &local[i][0], &remote[i][0]);
  }

  // Get the second share of each bit as computed by the other party
  exchange_shares_array_u(local_bits, remote_bits, numElements);

  for (int i=0; i<numElements; i++) {
    // Unpack the length/2 MSBs and store them at the last level
    unsigned long long tmp = ( remote_bits[i] >> (share_length/2) );
    memcpy(&remote[i][share_length-1], &tmp, sizeof(unsigned long long));
    // Unpack the rest and update odd levels
    for (int j=1; j<share_length-1; j+=2) {
      remote[i][j] = get_bit_u(remote_bits[i], j/2);
    }
  }

  /** REMAINING ROUNDS **/
  int rounds = (int) log2(share_length/2);
  char **levels = allocate_2D_byte_array(numElements, share_length/2);  // max 'length' levels per pair

  // Initialize level cache
  for (int i=0; i<numElements; i++) {
    for (int j=0; j<share_length/2; j++) {
      levels[i][j] = -1;
    }
  }

  int bits_left=share_length, bit_count;
  for (int r=1; r<=rounds; r++) {
    bits_left /= 2;
    for (int i=0; i<numElements; i++) {
      bit_count = 0;
      local_bits[i] = gr_round_c_char(r, bits_left, share_length,
                                 &local[i][0], &remote[i][0],
                                 &levels[i][0], &bit_count);
    }

    // Exchange all bits of the current round and unpack accordingly
    exchange_shares_array_u(local_bits, remote_bits, numElements);

    // Unpack bits of last level
    for (int i = 0; i < numElements; i++) {
      unsigned long long tmp = ( remote_bits[i] >> bit_count );
      memcpy(&remote[i][share_length-1], &tmp, sizeof(BShare));
      // Unpack the rest and reset level cache for next round
      int l=0;
      while ((levels[i][l] >= 0) & (l < share_length/2)) {
        remote[i][(int)levels[i][l]] = get_bit_u(remote_bits[i], l);
        levels[i][l++] = -1;  // Reset for next round
      }
    }
  }

  free(local_bits); free(remote_bits); free(levels);

  // One bitshare for each greater() comparison
  BitShare mask = 1;
  for (int i=0; i<numElements; i++) {
    result[i] = 0;
    // Do a final XOR of all levels
    for (int j=0; j<share_length-1; j++) {
      result[i] ^= local[i][j];
    }
    // XOR with last level
    result[i] ^= *((unsigned long long*) &local[i][share_length-1]);
    result[i] &= mask;
  }

  free(local); free(remote);

}

void select_eq_batch_const(BShareTable input, int leftcol, BShare cshare1,
                           BShare cshare2,  BShare result[]) {

  BShare** c = input.content;
  int len = input.numRows;
  BShare *res2 = (BShare *) malloc(len*sizeof(BShare)); // remote shares
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);

  // compute bitwise x^y^1
  for (long i=0; i<len; i++) {
   result[i] = c[i][leftcol] ^ cshare1 ^ (~(BShare)0); // local share;
   res2[i] = c[i][leftcol+1] ^ cshare2 ^ (~(BShare)0); // remote share
  }

  // The result is stored in the (numbits/2) rightmost bits of res, res2 elements
  for (int l=0; l<numlevels; l++) {

   for (long i=0; i<len; i++) {
     result[i] = eq_b_level2(numbits >> l, result[i], res2[i]);
   }

   // Exchange results of logical and, except for the final round
   if (l != numlevels-1) {
     exchange_shares_array(result, res2, len);
   }

  }
  // The local share of the final result is stored in res.
  // We need to call exchange again before using it in a subsequent operation.
  free(res2);
}

// Batched RCA with a const
void boolean_addition_batch_const(BShareTable input, int leftcol,
                                  BShare cshare1, BShare cshare2,
                                  BShare *res) {
  int numElements = input.numRows;
  BShare** t = input.content;
  BShare xor_b=0;
  BShare **c1 = allocate_2D_table(numElements, 3);
  BShare **c2 = allocate_2D_table(numElements, 3);

  BShare *carry1 = (BShare *) malloc(numElements * sizeof(BShare));
  BShare *carry2 = (BShare *) malloc(numElements * sizeof(BShare));
  BShare mask = 1;

  // Add bits ingoring carries
  for (int i=0; i<numElements; i++) {
    res[i] = t[i][leftcol]^cshare1;
    // Compute carry at position 0 (x_0 AND y_0)
    carry1[i] = and_b(get_bit(t[i][leftcol], 0), get_bit(t[i][leftcol+1], 0),
                   get_bit(cshare1, 0), get_bit(cshare2, 0),
                   get_next_rb()) & mask;
  }

  exchange_shares_array(carry1, carry2, numElements);
  for (int i=0; i<numElements; i++) {
    carry2[i] &= mask;
  }

  int bsize = sizeof(BShare)*8;

  // Compute carry at position j, ignoring the possible overflow carry
  for (int j=1; j<bsize - 1; j++) {
    for (int i=0; i<numElements; i++) {
      // x_j AND y_j
      c1[i][0] = and_b(get_bit(t[i][leftcol], j), get_bit(t[i][leftcol+1], j),
                    get_bit(cshare1, j), get_bit(cshare2, j),
                    get_next_rb());
      // x_j AND c_{j-1}
      c1[i][1] = and_b(get_bit(t[i][leftcol], j), get_bit(t[i][leftcol+1], j),
                    get_bit(carry1[i], j-1),
                    get_bit(carry2[i], j-1),
                    get_next_rb());
      // y_i AND c_{i-1}
      c1[i][2] =  and_b(get_bit(cshare1, j), get_bit(cshare2, j),
                     get_bit(carry1[i], j-1),
                     get_bit(carry2[i], j-1),
                     get_next_rb());
    }
    // Get carry shares from the other party
    // except for the last round
    if (j < bsize-2) {
      exchange_shares_array(&c1[0][0], &c2[0][0], numElements*3);
    }
    for (int i=0; i<numElements; i++) {
      // Store locally computed bit in carry1
      xor_b = (c1[i][0] ^ c1[i][1] ^ c1[i][2]);
      carry1[i] |= ((xor_b & mask) << j);
      // Store other party's bit share in carry2
      xor_b = (c2[i][0] ^ c2[i][1] ^ c2[i][2]);
      carry2[i] |= ((xor_b & mask) << j);
    }
  }

  // Add carries
  for (int i=0; i<numElements; i++) {
    res[i] ^= (carry1[i] << 1);
  }
  free(c1); free(c2); free(carry1); free(carry2);
}

// Group-by-sum that uses a ripple-carry adder
// The given BShareTable must be sorted on the group-by key(s)
void group_by_sum_rca(BShareTable* table, unsigned* key_indices,
                                          unsigned num_keys) {
  BShare** c = table->content;
  BShare max=0xFFFFFFFFFFFFFFFF;
  BShare rnums[187];  // Used by RCA

  // Scan table and update sums
  int len = table->numCols/2 + 1;  // Number of attributes plus new sum
  BShare local[len], remote[len];
  BShare bs1, bs2;
  BShare local_bits[num_keys], remote_bits[num_keys];
  // The indices of the shares for the sum attribute
  int cnt_idx_1 = table->numCols-2,   // First share
      cnt_idx_2 = table->numCols-1;   // Second share
  for (int i=0, k=0; i<table->numRows-1; i++, k+=2) {
    // Compute a single bit bs that denotes whether the adjacent rows c[i],
    // c[i+1] are in the same group
    for (int idx=0; idx<num_keys; idx++) {
      unsigned key_index = key_indices[idx];
      // Set local equality bit share
      local_bits[idx]  = eq_b(c[i][key_index], c[i][key_index+1],
                              c[i+1][key_index], c[i+1][key_index+1]);
    }
    // Fetch remote equality bit shares
    exchange_shares_array(local_bits, remote_bits, num_keys);   // 1 round
    // AND all equality bits to compute the final bit bs
    and_b_all(local_bits, remote_bits, num_keys);   // log(num_keys) rounds
    // bs denotes whether the adjacent rows c[i], c[i+1] are in the same group
    bs1 = local_bits[0];
    bs2 = remote_bits[0];
    BShare b1 = -bs1;  // Set all bits equal to LSB of bs1
    BShare b2 = -bs2;  // Set all bits equal to LSB of bs2
    // Compute new_s[i] = b * dummy_row + (1-b) * row1
    for (int j=0; j<table->numCols; j+=2) {
      local[j/2] = and_b(b1, b2, max, max, get_next_rb());
      local[j/2] ^= and_b(~b1, ~b2, c[i][j], c[i][j+1], get_next_rb());
    }
    // Compute new_sum = b*(sum[i] + sum[i+1]) + (1-b)*sum[i+1]
    BShare local_sum, remote_sum;
    get_next_rb_array(rnums, 187);
    boolean_addition(c[i][cnt_idx_1], c[i][cnt_idx_2],      // L rounds
                     c[i+1][cnt_idx_1], c[i+1][cnt_idx_2],
                     &local_sum, &remote_sum, rnums);
    local[len-1] = and_b(b1, b2, local_sum, remote_sum, get_next_rb());
    local[len-1] ^= and_b(~b1, ~b2, c[i+1][cnt_idx_1], c[i+1][cnt_idx_2],
                          get_next_rb());
    // Fetch remote boolean shares
    exchange_shares_array(local, remote, len);    // 1 round
    // Set c[i] = new_s[i]
    for (int j=0; j<table->numCols; j+=2) {
      c[i][j] = local[j/2];
      c[i][j+1] = remote[j/2];
    }
    // Set new_sum
    c[i+1][cnt_idx_1] = local[len-1];
    c[i+1][cnt_idx_2] = remote[len-1];
  }
}

// Group-by-sum that uses a ripple-carry adder and takes into account the given
// selected bits. The input BShareTable must be sorted on the group-by key(s)
void group_by_sum_rca_sel(BShareTable* table, BShare* selected,
                          unsigned* key_indices, unsigned num_keys) {
  BShare** c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  BShare rnums[187];  // Used by RCA

  // Fetch the second boolean share of each 'selected' bit -- 1 round
  AShare *remote_selected = (AShare *) malloc((table->numRows)*sizeof(BShare));
  assert(remote_selected!=NULL);
  exchange_shares_array(selected, remote_selected, table->numRows);

  // Scan table and update sums
  int len = table->numCols/2 + 3;  // Number of attributes plus the new sum
                                   // plus the old and new 'selected' bits
  BShare local[len], remote[len];
  BShare bs1, bs2, bn1, bn2, bk1, bk2, bh1, bh2, bj1, bj2;  // Bit shares
  BShare local_bits[num_keys], remote_bits[num_keys];
  // The indices of the shares for the sum attribute
  int cnt_idx_1 = table->numCols-2,   // First share
      cnt_idx_2 = table->numCols-1;   // Second share
  for (int i=0, k=0; i<table->numRows-1; i++, k+=2) {
    // Compute a single bit bs that denotes whether the adjacent rows c[i],
    // c[i+1] are in the same group
    for (int idx=0; idx<num_keys; idx++) {
      unsigned key_index = key_indices[idx];
      // Set local equality bit share
      local_bits[idx]  = eq_b(c[i][key_index], c[i][key_index+1],
                              c[i+1][key_index], c[i+1][key_index+1]);
    }
    // Fetch remote equality bit shares
    exchange_shares_array(local_bits, remote_bits, num_keys);   // 1 round
    // AND all equality bits to compute the final bit bs
    and_b_all(local_bits, remote_bits, num_keys);   // log(num_keys) rounds
    // bs denotes whether the adjacent rows c[i], c[i+1] are in the same group
    bs1 = local_bits[0];
    bs2 = remote_bits[0];
    // Compute bn = (bs OR NOT selected_b[i]) = NOT(NOT bs AND selected_b[i])
    // The composite bit bn is used in multiplexing the i-th tuple
    bn1 = bs1 ^ mask,
    bn2 = bs2 ^ mask;
    bn1 = and_b(bn1, bn2, selected[i] & mask, remote_selected[i] & mask,
                get_next_rb()) ^ mask;
    bn2 = exchange_shares(bn1);         // 1 round
    bn1 &= mask;                        // Keep LSB only
    bn2 &= mask;                        // Keep LSB only
    BShare b1 = -bn1;  // Set all bits equal to LSB of bs1
    BShare b2 = -bn2;  // Set all bits equal to LSB of bs2
    // Compute new_s[i] = b * dummy_row + (1-b) * row1
    for (int j=0; j<table->numCols; j+=2) {
      local[j/2] = and_b(b1, b2, max, max, get_next_rb());
      local[j/2] ^= and_b(~b1, ~b2, c[i][j], c[i][j+1], get_next_rb());
    }
    // Compute bk = bs AND selected[i]
    bk1 = and_b(bs1, bs2, selected[i], remote_selected[i],
                get_next_rb());
    bk2 = exchange_shares(bk1);  // 1 round
    bk1 &= mask;                        // Keep LSB only
    bk2 &= mask;                        // Keep LSB only
    // Compute bh = bs AND selected[i] AND selected[i+1]
    // ("same group and both selected")
    bh1 = and_b(bk1, bk2, selected[i+1], remote_selected[i+1],
                get_next_rb());
    bh2 = exchange_shares(bh1);  // 1 round
    bh1 &= mask;                        // Keep LSB only
    bh2 &= mask;                        // Keep LSB only
    bh1 = -bh1;                         // Set all bits equal to LSB
    bh2 = -bh2;                         // Set all bits equal to LSB
    // Compute bj = bs AND selected[i] AND (NOT selected[i+1])
    // ("same group and first selected")
    bj1 = and_b(bk1, bk2, selected[i+1] ^ mask,
                       remote_selected[i+1] ^ mask, get_next_rb());
    bj2 = exchange_shares(bj1);  // 1 round
    bj1 &= mask;                        // Keep LSB only
    bj2 &= mask;                        // Keep LSB only
    bj1 = -bj1;                         // Set all bits equal to LSB
    bj2 = -bj2;                         // Set all bits equal to LSB
    // Compute new_sum = bh*(sum[i] + sum[i+1]) + bj*sum[i] + (1-bs)*sum[i+1]
    BShare local_sum, remote_sum;
    get_next_rb_array(rnums, 187);
    boolean_addition(c[i][cnt_idx_1], c[i][cnt_idx_2],      // L rounds
                     c[i+1][cnt_idx_1], c[i+1][cnt_idx_2],
                     &local_sum, &remote_sum, rnums);
    local[len-3] = and_b(bh1, bh2, local_sum, remote_sum, get_next_rb());
    local[len-3] ^= and_b(bj1, bj2, c[i][cnt_idx_1], c[i][cnt_idx_2],
                          get_next_rb());
    local[len-3] ^= and_b(~(-bs1), ~(-bs2), c[i+1][cnt_idx_1],
                          c[i+1][cnt_idx_2], get_next_rb());
    // Compute bn * max + (1-bn)*selected[i] to allow masking old 'selected' bit
    local[len-2] = and_b(bn1, bn2, max, max, get_next_r());
    local[len-2] ^= and_b(~bn1, ~bn2, selected[i], remote_selected[i],
                        get_next_rb());
    // Compute cond = bs AND selected[i] to propagate 'selected' bit
    BShare cond1 = and_b(bs1, bs2, selected[i], remote_selected[i],
                         get_next_rb());
    BShare cond2 = exchange_shares(cond1);  // 1 round
    cond1 = -cond1;                         // Set all bits equal to LSB
    cond2 = -cond2;                         // Set all bits equal to LSB
    // Compute cond * selected[i] + (1-cond)*selected[i+1]
    local[len-1] = and_b(cond1, cond2, selected[i], remote_selected[i],
                         get_next_rb());
    local[len-1] ^= and_b(~cond1, ~cond2,
                          selected[i+1], remote_selected[i+1],
                          get_next_rb());
    // Fetch remote boolean shares
    exchange_shares_array(local, remote, len);    // 1 round
    // Set c[i] = new_s[i]
    for (int j=0; j<table->numCols; j+=2) {
      c[i][j] = local[j/2];
      c[i][j+1] = remote[j/2];
    }
    // Propagate 'selected' bit
    selected[i+1] = local[len-1];
    remote_selected[i+1] = remote[len-1];
    // Update i-th 'selected' bit
    selected[i] = local[len-2];
    remote_selected[i] = remote[len-2];
    // Set new_sum
    c[i+1][cnt_idx_1] = local[len-3];
    c[i+1][cnt_idx_2] = remote[len-3];
  }
  // Make sure we mask the last row if it's not selected
  // 1. Compute composite bit bl = NOT(bs) AND NOT(selected)
  BShare bl1 = and_b(bs1 ^ mask, bs2 ^ mask,
                     selected[table->numRows-1] ^ mask,
                     remote_selected[table->numRows-1] ^ mask,
                     get_next_rb())
                   & mask;
  BShare bl2 = exchange_shares(bl1);            // 1 round
  BShare b1 = -bl1;
  BShare b2 = -bl2;
  // 2. Compute new_c[i] = bl * dummy_row + (1-bl) * last_row
  for (int j=0; j<table->numCols; j+=2) {
    local[j/2] = and_b(b1, b2, max, max, get_next_rb());
    local[j/2] ^= and_b(~b1, ~b2, c[table->numRows-1][j],
                        c[table->numRows-1][j+1], get_next_rb());
  }
  exchange_shares_array(local, remote, len);    // 1 round
  // 3. Multiplex
  for (int j=0; j<table->numCols; j+=2) {
    c[table->numRows-1][j] = local[j/2];
    c[table->numRows-1][j+1] = remote[j/2];
  }
  free(remote_selected);
  // Final shuffle iff not followed by ORDER_BY
}

void group_by_sum_rca_sel_odd_even(BShareTable* table, int batch_size, 
                                   BShare* selected_b,
                                   unsigned* key_indices, unsigned num_keys) {
  // Make sure the number of rows is a power of two
  assert(ceil(log2(table->numRows)) == floor(log2(table->numRows)));
  // Make sure batch_size is at most equal to the input size
  batch_size = batch_size > table->numRows ? table->numRows : batch_size;
  BShare** c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Fetch the second boolean share of each 'selected' bit (1 round)
  BShare *remote_selected_b = (BShare *) malloc((table->numRows)*sizeof(BShare));
  assert(remote_selected_b!=NULL);
  exchange_shares_array(selected_b, remote_selected_b, table->numRows);
  // Allocate memory for local and remote shares
  BShare  *local_bits= (BShare *) malloc(batch_size*num_keys*sizeof(BShare)),  // local shares of equality bits
          *remote_bits= (BShare *) malloc(batch_size*num_keys*sizeof(BShare)); // remote shares of equality bits
  assert(local_bits!=NULL); assert(remote_bits!=NULL);
  // Allocate memory for local and remote shares
  BShare  *local_bns= (BShare *) malloc(batch_size*num_keys*sizeof(BShare)),  // local shares of bn bits
          *remote_bns= (BShare *) malloc(batch_size*num_keys*sizeof(BShare)); // remote shares of bn bits
  assert(local_bns!=NULL); assert(remote_bns!=NULL);
  BShare  *local_sums= (BShare *) malloc(batch_size*sizeof(BShare)),  // local shares of sums
          *remote_sums= (BShare *) malloc(batch_size*sizeof(BShare)); // remote shares of sums
  assert(local_sums!=NULL); assert(remote_sums!=NULL); 
  int len = table->numCols/2 + 1;  // Number of attributes plus new count
  BShare  *local= (BShare *) malloc(batch_size*len*sizeof(BShare)),  // local shares of multiplexed rows
          *remote= (BShare *) malloc(batch_size*len*sizeof(BShare)); // remote shares of multiplexed rows
  assert(local!=NULL); assert(remote!=NULL);
  // The indices of the shares for the sum attribute
  int cnt_idx_1 = table->numCols-2;   // First share index
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int rounds = log2(table->numRows);  // Number of rounds for odd-even aggregation
  int dist = table->numRows;  // Distance between elements to compare at each round
  int num_pairs, start, step, end, offset, total_comparisons, num_comparisons=0;
  // Use an odd-even aggregation and count rows by adding the 'selected' bits
  for (int r=0; r<rounds; r++){
    dist /= 2;
    total_comparisons = table->numRows-dist;
    start = 0;
    end = (batch_size <= total_comparisons ? batch_size : total_comparisons);
    num_pairs = end - start;
    do {  // For all batches 
      // For all pairs of elements in the current batch
      for (int i=start, k=0; i<end; i++, k+=num_keys) { 
        for (int idx=0; idx<num_keys; idx++) { // For all group-by keys
            unsigned key_index = key_indices[idx];
            // Compute bitwise x^y^1
            local_bits[k+idx] = c[i][key_index] ^ c[i+dist][key_index] ^ (~(BShare)0); 
            remote_bits[k+idx] = c[i][key_index+1] ^ c[i+dist][key_index+1] ^ (~(BShare)0); 
        }
      }
      // Compute equality bits in bulk
      eq_bulk(numlevels, numbits, local_bits, remote_bits, num_pairs*num_keys); // numlevels-1 rounds
      exchange_shares_array(local_bits, remote_bits, num_pairs*num_keys);  // 1 round
      // AND all equality bits per pair to compute the final bit bs
      and_b_all_group(local_bits, remote_bits, num_pairs, num_keys);  // log(num_keys) rounds
      // Compute masking bit bn = (bs OR NOT selected_b[i+dist]) = NOT(NOT bs AND selected_b[i+dist])
      for (int i=start, k=0; i<end; i++, k++) {
        local_bns[k] = local_bits[k] ^ mask;
        remote_bns[k] = remote_bits[k] ^ mask; 
        local_bns[k] = and_b(local_bns[k], remote_bns[k], selected_b[i+dist] & mask, remote_selected_b[i+dist] & mask,
                              get_next_rb()) ^ mask;
      }
      exchange_shares_array(local_bns, remote_bns, num_pairs);  // 1 round
      // Compute new_c[i+dist] = bn * dummy_row + (1-bn) * c[i+dist] to multiplex rows
      for (int i=start, k=0; i<end; i++, k++) {
        local_bns[k] &= mask;  // Keep LSB only
        local_bns[k] = -local_bns[k];  // Set all bits equal to the LSB
        remote_bns[k] &= mask;  // Keep LSB only
        remote_bns[k] = -remote_bns[k];  // Set all bits equal to the LSB
        offset = k*len;
        for (int j=0; j<table->numCols; j+=2) {
          local[offset+j/2] = and_b(local_bns[k], remote_bns[k], max, max, get_next_rb());
          local[offset+j/2] ^= and_b(~local_bns[k], ~remote_bns[k], c[i+dist][j], c[i+dist][j+1], get_next_rb());
        }
      }
      // Compute sum = c[i][cnt_idx] + c[i+dist][cnt_idx]
      boolean_addition_batch2(c, local_sums, start, dist, num_pairs, cnt_idx_1);  // logL rounds
      exchange_shares_array(local_sums, remote_sums, num_pairs);  // 1 round
      // Compute bu = bs AND selected[i] AND selected[i+dist]
      for (int i=start, k=0; i<end; i++, k++) {
        local_bits[k] = and_b(local_bits[k], remote_bits[k], 
                        selected_b[i], remote_selected_b[i],
                        get_next_rb());
      }
      exchange_shares_array(local_bits, remote_bits, num_pairs);  // 1 round
      for (int i=start, k=0; i<end; i++, k++) {
        local_bits[k] = and_b(local_bits[k], remote_bits[k], 
                        selected_b[i+dist], remote_selected_b[i+dist],
                        get_next_rb());
      }
      exchange_shares_array(local_bits, remote_bits, num_pairs);  // 1 round

      for (int i=start, k=0; i<end; i++, k++) {
        local_bits[k] &= mask;
        remote_bits[k] &= mask;
        local_bits[k] = -local_bits[k];
        remote_bits[k] = -remote_bits[k];
        // Compute new_cnt = bu*sum + (1-bu)*c[i][cnt_idx]
        offset = (k+1) * len;
        local[offset-1] = and_b(local_bits[k], remote_bits[k], local_sums[k], remote_sums[k],
                                get_next_rb());
        local[offset-1] ^= and_b(~local_bits[k], ~remote_bits[k], c[i][cnt_idx_1],
                                 c[i][cnt_idx_1+1], get_next_rb());
      }
      // Fetch remote boolean and arithmetic shares in a single round
      // NOTE (john): This works because BShare and AShare are both of the same type
      exchange_shares_array(local, remote, num_pairs*len);    // 1 round
      // Update rows in place
      for (int i=start, k=0; i<end; i++, k++) {
        offset = k*len;
        // Set c[i+dist] = new_c[i+dist]
        for (int j=0; j<table->numCols; j+=2) {
          c[i+dist][j] = local[offset+j/2];
          c[i+dist][j+1] = remote[offset+j/2];
        }
        // Set new_cnt
        c[i][cnt_idx_1] = local[offset+len-1];
        c[i][cnt_idx_1+1] = remote[offset+len-1];
      }
      num_comparisons += num_pairs;
      // Update batch boundaries
      start = end;
      step = end + batch_size;
      end = (step <= total_comparisons ? step : total_comparisons);
      num_pairs = end - start;
    } while (start<end);
  }
  // Do a final pass to do the necessary multiplexing  
  // Compute bit = selected[0]==0 to multiplex first row
  BShare bit, rbit;
  local[0] = selected_b[0] ^ 0 ^ (~(BShare)0); 
  remote[0] = remote_selected_b[0] ^ 0 ^ (~(BShare)0);
  for (int start=0; start<table->numRows; start+=batch_size) {
    int end = (start+batch_size <= table->numRows ? start+batch_size : table->numRows);
    int comp = end - start; 
    offset = start==0 ? 1 : 0;
    // Compute b = c[start][0]==max to multiplex the counts
    for (int i=start+offset, k=offset; i<end; i++, k++) {
      local[k] = c[i][0] ^ max ^ (~(BShare)0); 
      remote[k] = c[i][1] ^ max ^ (~(BShare)0);
    }
    // Compute equality bits in bulk
    eq_bulk(numlevels, numbits, local, remote, comp); // numlevels-1 rounds
    exchange_shares_array(local, remote, comp);  // 1 round
    if (start==0) {
      bit = local[0] & mask;  // Keep LSB only
      bit = -bit;  // Set all bits equal to the LSB
      rbit = remote[0] & mask;  // Keep LSB only
      rbit = -rbit;  // Set all bits equal to the LSB
      // Compute c[0] = bit*max + (1-bit)*c[0]
      for (int j=0; j<table->numCols; j+=2) {
        local[j/2] = and_b(bit, rbit, max, max, get_next_rb());
        local[j/2] ^= and_b(~bit, ~rbit, c[0][j], c[0][j+1], get_next_rb());
      }
      offset = table->numCols/2;
    }
    for (int i=start, k=0; i<end; i++, k++) {
      bit = local[k] & mask;  // Keep LSB only
      bit = -bit;  // Set all bits equal to the LSB
      rbit = remote[k] & mask;  // Keep LSB only
      rbit = -rbit;  // Set all bits equal to the LSB
      // Compute selected[i] = b*max + (1-b)*selected[i]
      local[k+offset] = and_b(bit, rbit, max, max, get_next_rb());
      local[k+offset] ^= and_b(~bit, ~rbit,
                               selected_b[i], remote_selected_b[i], get_next_rb());
    }
    exchange_shares_array(local, remote, comp+offset); // 1 round
    if (start==0) {  // Update first row
      for (int j=0; j<table->numCols; j+=2) {
        c[0][j] = local[j/2];
        c[0][j+1] = remote[j/2];
      }
    }
    // Update counts
    for (int i=start, k=0; i<end; i++, k++) {
      selected_b[k] = local[k+offset];
      remote_selected_b[k] = remote[k+offset];
    }
  }
  // TODO (john): Final shuffle iff not followed by ORDER_BY
  free(remote_selected_b);
  free(local_sums); free(remote_sums);
  free(local_bns); free(remote_bns);
  free(local_bits); free(remote_bits);
  free(local); free(remote);                             
}

void group_by_sum_rca_odd_even(BShareTable* table, int batch_size, 
                               unsigned* key_indices, unsigned num_keys) {
  // Make sure the number of rows is a power of two
  assert(ceil(log2(table->numRows)) == floor(log2(table->numRows)));
  // Make sure batch_size is at most equal to the input size
  batch_size = batch_size > table->numRows ? table->numRows : batch_size;
  BShare** c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Allocate memory for local and remote shares
  BShare  *local_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)),  // local shares of equality bits
          *remote_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)); // remote shares of equality bits
  assert(local_bits!=NULL); assert(remote_bits!=NULL);
  BShare  *local_sums= (BShare  *) malloc(batch_size*sizeof(BShare)),  // local shares of sums
          *remote_sums= (BShare  *) malloc(batch_size*sizeof(BShare)); // remote shares of sums
  assert(local_sums!=NULL); assert(remote_sums!=NULL); 
  int len = table->numCols/2 + 1;  // Number of attributes plus new count
  BShare  *local= (BShare  *) malloc(batch_size*len*sizeof(BShare)),  // local shares of multiplexed rows
          *remote= (BShare  *) malloc(batch_size*len*sizeof(BShare)); // remote shares of multiplexed rows
  assert(local!=NULL); assert(remote!=NULL);
  // The indices of the shares for the sum attribute
  int cnt_idx_1 = table->numCols-2;   // First share index
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int rounds = log2(table->numRows);  // Number of rounds for odd-even aggregation
  int dist = table->numRows;  // Distance between elements to compare at each round
  int num_pairs, start, step, end, offset, total_comparisons, num_comparisons=0;
  // Use an odd-even aggregation and count rows by adding the 'selected' bits
  for (int r=0; r<rounds; r++){
    dist /= 2;
    total_comparisons = table->numRows-dist;
    start = 0;
    end = (batch_size <= total_comparisons ? batch_size : total_comparisons);
    num_pairs = end - start;
    do {  // For all batches 
      // For all pairs of elements in the current batch
      for (int i=start, k=0; i<end; i++, k+=num_keys) { 
        for (int idx=0; idx<num_keys; idx++) { // For all group-by keys
            unsigned key_index = key_indices[idx];
            // Compute bitwise x^y^1
            local_bits[k+idx] = c[i][key_index] ^ c[i+dist][key_index] ^ (~(BShare)0); 
            remote_bits[k+idx] = c[i][key_index+1] ^ c[i+dist][key_index+1] ^ (~(BShare)0); 
        }
      }
      // Compute equality bits in bulk
      eq_bulk(numlevels, numbits, local_bits, remote_bits, num_pairs*num_keys); // numlevels-1 rounds
      exchange_shares_array(local_bits, remote_bits, num_pairs*num_keys);  // 1 round
      // AND all equality bits per pair to compute the final bit bs
      and_b_all_group(local_bits, remote_bits, num_pairs, num_keys);  // log(num_keys) rounds
      // Compute new_c[i+dist] = bs*dummy_row + (1-bs)*c[i+dist] to multiplex rows
      for (int i=start, k=0; i<end; i++, k++) {
        local_bits[k] &= mask;  // Keep LSB only
        local_bits[k] = -local_bits[k];  // Set all bits equal to the LSB
        remote_bits[k] &= mask;  // Keep LSB only
        remote_bits[k] = -remote_bits[k];  // Set all bits equal to the LSB
        offset = k*len;
        for (int j=0; j<table->numCols; j+=2) {
          local[offset+j/2] = and_b(local_bits[k], remote_bits[k], max, max, get_next_rb());
          local[offset+j/2] ^= and_b(~local_bits[k], ~remote_bits[k], c[i+dist][j], c[i+dist][j+1], get_next_rb());
        }
      }
      // Compute sum = c[i][cnt_idx] + c[i+dist][cnt_idx] 
      boolean_addition_batch2(c, local_sums, start, dist, num_pairs, cnt_idx_1);
      exchange_shares_array(local_sums, remote_sums, num_pairs);
      // Compute new_cnt = bs*sum + (1-bs)*c[i][cnt_idx] 
      for (int i=start, k=0; i<end; i++, k++) {
        offset = (k+1) * len;
        local[offset-1] = and_b(local_bits[k], remote_bits[k], local_sums[k], remote_sums[k],
                              get_next_rb());
        local[offset-1] ^= and_b(~local_bits[k], ~remote_bits[k], c[i][cnt_idx_1],
                                 c[i][cnt_idx_1+1], get_next_rb());
      }
      // Fetch remote boolean and arithmetic shares in a single round
      // NOTE (john): This works because BShare and AShare are both of the same type
      exchange_shares_array(local, remote, num_pairs*len);    // 1 round
      // Update rows in place
      for (int i=start, k=0; i<end; i++, k++) {
        offset = k*len;
        // Set c[i+dist] = new_c[i+dist]
        for (int j=0; j<table->numCols; j+=2) {
          c[i+dist][j] = local[offset+j/2];
          c[i+dist][j+1] = remote[offset+j/2];
        }
        // Set c[i][cnt_idx] = new_cnt
        c[i][cnt_idx_1] = local[offset+len-1];
        c[i][cnt_idx_1+1] = remote[offset+len-1];
      }
      num_comparisons += num_pairs;
      // Update batch boundaries
      start = end;
      step = end + batch_size;
      end = (step <= total_comparisons ? step : total_comparisons);
      num_pairs = end - start;
    } while (start<end);
  }
  // Do a final pass to properly mask rows
  for (int start=0; start<table->numRows; start+=batch_size) {
    int end = (start+batch_size <= table->numRows ? start+batch_size : table->numRows);
    int comp = end - start;
    // Compute b = c[start][0]==max that denotes whether the row has been masked (1) or not (0)
    for (int i=start, k=0; i<end; i++, k++) {
      local_bits[k] = c[i][0] ^ max ^ (~(BShare)0); 
      remote_bits[k] = c[i][1] ^ max ^ (~(BShare)0);
    }
    // Compute equality bits in bulk
    eq_bulk(numlevels, numbits, local_bits, remote_bits, comp); // numlevels-1 rounds
    exchange_shares_array(local_bits, remote_bits, comp);  // 1 round
    // Compute c[i] = b*max + (1-b)*c[i]
    for (int i=start, k=0; i<end; i++, k++) {
      local_bits[k] &= mask;
      remote_bits[k] &= mask;
      local_bits[k] = -local_bits[k];
      remote_bits[k] = -remote_bits[k];
      offset = k * len;
      for (int j=0; j<table->numCols; j+=2) {
        local[offset+j/2] = and_b(local_bits[k], remote_bits[k], max, max, get_next_rb());
        local[offset+j/2] ^= and_b(~local_bits[k], ~remote_bits[k],
                                    c[i][j], c[i][j+1], get_next_rb());
      }
    }
    exchange_shares_array(local, remote, comp*len); // 1 round
    // Update rows
    for (int i=start, k=0; i<end; i++, k++) {
      offset = k * len;
      for (int j=0; j<table->numCols; j+=2) {
        c[i][j] = local[offset+j/2];
        c[i][j+1] = remote[offset+j/2];
      }      
    }
  }
  // TODO (john): Final shuffle iff not followed by ORDER_BY
  free(local_bits); free(remote_bits);
  free(local_sums); free(remote_sums);
  free(local); free(remote);                        
}


PRIVATE unsigned long long gr_round_a(BShare x1, BShare x2, BShare y1, BShare y2, int length) {
  // Compute (x_i ^ y_i)
  BShare xor1 = x1 ^ y1;
  BShare xor2 = x2 ^ y2;

  unsigned long long last_and = 0, last_ands = 0;
  const BShare mask=1;
  int index;  // The bit index (index=0 for the LSB)

  for (int i=0; i<length-1; i++) {
    // Compute ((x_{length-i-1} ^ y_{length-i-1}) AND x_{length-i-1})
    index = length-i-1;
    last_and = and_b(get_bit(xor1, index), get_bit(xor2, index),
                     get_bit(x1, index), get_bit(x2, index), get_next_rb())
                    & mask;
    // Pack result bit in last_ands to send all together
    last_ands |= ( last_and << index );
  }

  // Set LSB to (x_0 AND ~y_0)
  last_ands |= and_b(get_bit(x1, 0), get_bit(x2, 0),
                          get_bit(y1, 0) ^ mask, get_bit(y2, 0) ^ mask,
                          get_next_rb()) & mask;

  return last_ands;
}

PRIVATE unsigned long long geq_round_a(BShare x1, BShare x2, BShare y1,
                                       BShare y2, int length) {
  // Compute (x_i ^ y_i)
  BShare xor1 = x1 ^ y1;
  BShare xor2 = x2 ^ y2;

  unsigned long long last_and = 0, last_ands = 0;
  const BShare mask=1;
  int index;  // The bit index (index=0 for the LSB)

  for (int i=0; i<length-1; i++) {
    // Compute ((x_{length-i-1} ^ y_{length-i-1}) AND x_{length-i-1})
    index = length-i-1;
    last_and = and_b(get_bit(xor1, index), get_bit(xor2, index),
                     get_bit(x1, index), get_bit(x2, index), get_next_rb())
                    & mask;
    // Pack result bit in last_ands to send all together
    last_ands |= ( last_and << index );
  }

  // Store ~(~x_0 AND y_0) at the last level (length-1)
  last_ands |= ( and_b(get_bit(x1, 0) ^ mask, get_bit(x2, 0) ^ mask,
                       get_bit(y1, 0), get_bit(y2, 0),
                       get_next_rb())
                      & mask ) ^ mask;

  // Set LSB to (x_0 AND ~y_0)
  // last_ands |= and_b(get_bit(x1, 0), get_bit(x2, 0),
  //                         get_bit(y1, 0) ^ mask, get_bit(y2, 0) ^ mask,
  //                         get_next_rb()) & mask;

  return last_ands;
}

// B. Compute next to last AND at odd levels as well as 1st round of pairwise
// ANDs at the last level. This step performs 'length' logical ANDs in total.
PRIVATE unsigned long long gr_round_b(BShare x1, BShare x2, BShare y1, BShare y2,
                                      int length, char local[], char remote[]) {
  // Compute ~(x_i ^ y_i)
  BShare not_xor1 = ~(x1^y1);
  BShare not_xor2 = ~(x2^y2);
  int index;
  const BShare mask=1;

  unsigned long long local_bits = 0;
  for (int i=1, j=0; i<length-1; i+=2, j++) { // For all odd levels (length/2)
    index = length-i;
    // Set ~(x_{length-i} ^ y_{length-i}) next to the last bit
    local[i] |= ( get_bit(not_xor1, index) << 1 );
    // Set ~(x_{length-i} ^ y_{length-i}) next to the last remote bit
    remote[i] |= ( get_bit(not_xor2, index) << 1 );
    // Compute next to last logical AND for level i
    local[i] = ( and_b(get_bit_u8(local[i], 0), get_bit_u8(remote[i], 0),
                       get_bit_u8(local[i], 1), get_bit_u8(remote[i], 1),
                       get_next_rb()) & mask );
    // Pack result bit in local_bits to send all together
    local_bits |= (local[i] << j);
  }
  // Compute first round of pairwise logical ANDs at the last level
  unsigned long long tmp = eq_b_level2(length,
        unset_lsbs(not_xor1, 1) | *((unsigned long long*) &local[length-1]),
        unset_lsbs(not_xor2, 1) | *((unsigned long long*)  &remote[length-1]));
  memcpy(&local[length-1], &tmp, sizeof(unsigned long long));

  // Pack the length/2 result bits in the vacant MSBs of local_bits
  // local_bits |= ( ((unsigned long long) local[length-1]) << (length/2) );
  local_bits |= ( (*((unsigned long long*) &local[length-1])) << (length/2) );

  return local_bits;
}

// C. Continue in a loop. Each round breaks into the the following steps:
//    1. Project every other bit of the last level to 2^r levels,
//       starting at level 2^r * (bits_left - p), where p is the bit's index
//       (p=0 for the LSB). The bit is copied next to the LSB of the
//       corresponding level.
//    2. Evaluate a logical AND between the projected bit and the LSB at
//       the corresponding level.
//    3. Evaluate the next round of pairwise ANDs at the last level.
PRIVATE unsigned long long gr_round_c(int i, int bits_left, int length,
                                      char local[], char remote[],
                                      int levels[], int *bit_count) {

  int current_level, num_levels;
  const BShare mask=1;
  unsigned long long to_send = 0;
  num_levels = (1 << i);
  // Project common bits to avoid redundant computation (and communication)
  for (int p=bits_left-1; p>0; p-=2) {
    current_level = num_levels * (bits_left - p);
    for (int j=0; j<num_levels; j++) {
      // Project bits from last level
      BShare l_tmp = *((unsigned long long*) &local[length-1]);
      BShare r_tmp = *((unsigned long long*) &remote[length-1]);
      local[current_level] |= ( get_bit(l_tmp, p) << 1 );
      remote[current_level] |= ( get_bit(r_tmp, p) << 1 );
      // Do the logical AND
      local[current_level] = and_b(get_bit_u8(local[current_level], 0),
                                   get_bit_u8(remote[current_level], 0),
                                   get_bit_u8(local[current_level], 1),
                                   get_bit_u8(remote[current_level], 1),
                                   get_next_rb()) & mask;
      // Pack the result
      to_send |= ( local[current_level] << (*bit_count) );
      // Cache level to unpack remote bit later
      levels[*bit_count] = current_level;
      (*bit_count)++;
      current_level++;
      if ( current_level == (length-1) ) break;
    }
  }
  // Process last level
  BShare tmp = eq_b_level2(bits_left, *((unsigned long long*) &local[length-1]),
                                    *((unsigned long long*) &remote[length-1]));
  memcpy(&local[length-1], &tmp, sizeof(unsigned long long));
  // Pack bits of the last level
  to_send |= ( *((unsigned long long*) &local[length-1]) << (*bit_count) );
  return to_send;
}

PRIVATE unsigned long long gr_round_c_char(int i, int bits_left, int length, char local[], char remote[],
                                      char levels[], int *bit_count) {

  int current_level, num_levels;
  const BShare mask=1;
  unsigned long long to_send = 0;
  num_levels = (1 << i);
  // Project common bits to avoid redundant computation (and communication)
  for (int p=bits_left-1; p>0; p-=2) {
    current_level = num_levels * (bits_left - p);
    for (int j=0; j<num_levels; j++) {
      // Project bits from last level
      BShare l_tmp = *((unsigned long long*) &local[length-1]);
      BShare r_tmp = *((unsigned long long*) &remote[length-1]);
      local[current_level] |= ( get_bit(l_tmp, p) << 1 );
      remote[current_level] |= ( get_bit(r_tmp, p) << 1 );
      // Do the logical AND
      local[current_level] = and_b(get_bit_u8(local[current_level], 0),
                                   get_bit_u8(remote[current_level], 0),
                                   get_bit_u8(local[current_level], 1),
                                   get_bit_u8(remote[current_level], 1),
                                   get_next_rb()) & mask;
      // Pack the result
      to_send |= ( local[current_level] << (*bit_count) );
      // Cache level to unpack remote bit later
      levels[*bit_count] = current_level;
      (*bit_count)++;
      current_level++;
      if ( current_level == (length-1) ) break;
    }
  }
  // Process last level
  BShare tmp = eq_b_level2(bits_left, *((unsigned long long*) &local[length-1]),
                                    *((unsigned long long*) &remote[length-1]));
  memcpy(&local[length-1], &tmp, sizeof(unsigned long long));
  // Pack bits of the last level
  to_send |= ( *((unsigned long long*) &local[length-1]) << (*bit_count) );
  return to_send;
}

// Used by group_by_sum_rca_odd_even() and group_by_sum_rca_sel_odd_even()
PRIVATE void eq_bulk(int num_levels, int num_bits, BShare* local,
                     BShare* remote, int length) {
  // For all rounds of the equality check (6 in total)
  for (int l=0; l<num_levels; l++) {
    // Apply l-th round for all equalities in batch
    for (int k=0; k<length; k++) {
      local[k] = eq_b_level2(num_bits >> l, local[k], remote[k]);
    }
    // Exchange results of logical and, except for the final round
    if (l != num_levels-1) {
      exchange_shares_array(local, remote, length);
    }
  }
}
