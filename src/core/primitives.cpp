#include "../../include/core/primitives.h"
#include "mpi.h"

#define PRIVATE static
#define XCHANGE_MSG_TAG 7

PRIVATE void compute_composite_2(BitShare**, BitShare**, BitShare**, BitShare**, int, int);
PRIVATE void compute_composite_3(BitShare**, BitShare**, BitShare**, BitShare**, int, int);
PRIVATE int get_next_index(int, int, int);
PRIVATE void shift_greater(const BShare, const BShare, int, int, BShare*, BShare*);

inline PRIVATE void shift_greater(const BShare x1, const BShare x2, int r, int len,
                                  BShare *shifted_1, BShare *shifted_2) {
    int part_size = ((int) 1 << r);
    int offset = (part_size >> 1), start;
    BShare b1, b2;
    *shifted_1 = x1;
    *shifted_2 = x2;
    for (int i=0; i<len; i+=part_size) {
        start = i + offset;
        b1 = get_bit(x1, start);
        b2 = get_bit(x2, start);
        for (int j=start-1; j>=i; j--) {
            *shifted_1 = unset_bit(*shifted_1,j) | (b1 << j);
            *shifted_2 = unset_bit(*shifted_2,j) | (b2 << j);
        }
    }
}

// Addition
AShare add(AShare x, AShare y) {
  return x + y;
}

// Multiplication (assumes three parties)
AShare mul(AShare x1, AShare x2, AShare y1, AShare y2, AShare rnum) {
  AShare z = (x1 * y1) + (x1 * y2) + (x2 * y1);
  // Add random number
  z += rnum;
  return z;
}

// Equality (assumes three parties)
AShare eq(AShare x1, AShare x2, AShare y1, AShare y2, AShare w1, AShare w2,
         AShare rnum) {
  AShare s1 = x1 - y1;
  AShare s2 = x2 - y2;
  // Compute (x-y)w
  return mul(s1, s2, w1, w2, rnum);
}

// Logical OR (assumes three parties)
AShare or_a(AShare x1, AShare x2, AShare y1, AShare y2, AShare rnum) {
  return mul(x1, x2, y1, y2, rnum);
}

// Logical AND (assumes three parties)
AShare and_a(AShare x1, AShare x2, AShare y1, AShare y2, AShare rnum1,
             AShare rnum2) {
  AShare x_sq = mul(x1, x2, x1, x2, rnum1);
  AShare y_sq = mul(y1, y2, y1, y2, rnum2);
  return add(x_sq, y_sq);
}

// Complement
BShare not_b(BShare x) {
  return ~x;
}

// Logical AND using boolean shares (assumes three parties)
BShare and_b(const BShare x1, const BShare x2, const BShare y1,
                    const BShare y2, const BShare rnum) {
  BShare z = (x1 & y1) ^ (x1 & y2) ^ (x2 & y1);
  // XOR random number
  z ^= rnum;
  return z;
}

// Logical AND using boolean shares for arrays
void and_b_array(const BShare *x1, const BShare *x2,
                        const BShare *y1, const BShare *y2,
                        const BShare *rnum, int len, BShare *res) {

  for (int i=0; i<len; i++) {
    res[i] = (x1[i] & y1[i]) ^ (x1[i] & y2[i]) ^ (x2[i] & y1[i]);
    res[i] ^= rnum[i];
  }
}

// Logical AND of all elements in the given array
void and_b_all(BShare *x1, BShare *x2, int len) {
  assert(len>0);
  while (len>1) {
    // For each pair of elements in the given array
    for (int i=0; i<len-1; i+=2) {
      x1[i/2] = (x1[i] & x1[i+1]) ^ (x1[i] & x2[i+1]) ^ (x2[i] & x1[i+1]);
      x1[i/2] ^= get_next_rb();
    }
    // Get remote shares
    exchange_shares_array(x1, x2, len/2);
    int m = len % 2;
    if (m != 0) {  // Shift last odd element
      x1[len/2] = x1[len-1];
      x2[len/2] = x2[len-1];
    }
    // Number of remaining elements
    len = len/2 + m;
  }
  // Keep LSB
  x1[0] &= (BShare) 1;
  x2[0] &= (BShare) 1;
}

// Logical AND of all elements per group in the given array
void and_b_all_group(BShare *x1, BShare *x2, int num_groups, int group_size) {
  assert(num_groups>0 && group_size>0);
  while (group_size>1) { // For each round
    int end = num_groups*group_size;  // Total number of elements in the current round 
    int offset = group_size/2;        // Number of freed slots per group in the current round
    int m = group_size % 2;
    for (int start=0, k=0; start<end; start+=group_size, k++) {  // For each group
      int idx = start - k*offset;  // The new group start index 
      // For each pair of elements in the group
      for (int i=0; i<group_size-1; i+=2) {
        x1[idx+i/2] = (x1[start+i] & x1[start+i+1]) ^ (x1[start+i] & x2[start+i+1]) ^ (x2[start+i] & x1[start+i+1]);
        x1[idx+i/2] ^= get_next_rb();
      }
      if (m != 0) {  // Shift last odd element
        x1[idx+group_size/2] = x1[start+group_size-1];
        x2[idx+group_size/2] = x2[start+group_size-1];
      }
    }
    // Number of remaining elements per group
    group_size = group_size/2 + m;
    // Get remote shares
    exchange_shares_array(x1, x2, group_size*num_groups); 
  }
  // Keep LSB
  for (int i=0; i<num_groups; i++) {
    x1[i] &= (BShare) 1;
    x2[i] &= (BShare) 1;
  }
}

// Logical AND using single-bit boolean shares (assumes three parties)
inline BitShare and_bit(BitShare x1, BitShare x2, BitShare y1, BitShare y2,
                        BitShare rnum) {
  BitShare z = (x1 & y1) ^ (x1 & y2) ^ (x2 & y1);
  // XOR random bit
  z ^= rnum;
  return z;
}

// array-based and_bit
void and_bit_array(const BitShare *x1, const BitShare *x2,
                        const BitShare *y1, const BitShare *y2,
                        const BitShare *rnum, int len, BitShare *res) {

  for (int i=0; i<len; i++) {
    res[i] = (x1[i] & y1[i]) ^ (x1[i] & y2[i]) ^ (x2[i] & y1[i]);
    res[i] ^= rnum[i];
  }
}

// Logical XOR
BShare xor_b(BShare x, BShare y) {
  return x ^ y;
}

// Boolean equality
BShare eq_b(BShare x1, BShare x2, BShare y1, BShare y2) {
    // compute bitwise x^y^1
    BShare lbs = x1 ^ y1 ^ ~((BShare) 0); // local share
    BShare rbs = x2 ^ y2 ^ ~((BShare) 0); // remote share
    int numbits = sizeof(BShare) * 8;
    int rounds = (int) log2(numbits); // Number of rounds for oblivious equality
    BShare shifted_lbs, shifted_rbs;
    for (int l=0; l<rounds-1; l++) {
        // Compute first round
        shift_greater(lbs, rbs, l, numbits, &shifted_lbs, &shifted_rbs);
        lbs = and_b(lbs, rbs, shifted_lbs, shifted_rbs, get_next_rb());
        rbs = exchange_shares(lbs);
    }
    // Do final round without exchanging shares
    shift_greater(lbs, rbs, rounds-1, numbits, &shifted_lbs, &shifted_rbs);
    lbs = and_b(lbs, rbs, shifted_lbs, shifted_rbs, get_next_rb());
    // Return one share of the final result (LSB).
    // Note: We need to call exchange again before using it in a subsequent operation.
    return lbs & (BShare) 1;
}

// Blocking boolean equality (with eager rnum generation)
BShare eq_b_sync(BShare x1, BShare x2, BShare y1, BShare y2) {

  // compute bitwise x^y^1
  BShare res1 = x1 ^ y1 ^ (~(BShare)0); // local share
  BShare res2 = x2 ^ y2 ^ (~(BShare)0); // remote share

  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);

  // The result is stored in the (numbits/2) rightmost bits of res1, res2
  for (int l=0; l<numlevels; l++) {
    res1 = eq_b_level2(numbits >> l, res1, res2);

    // Exchange results of logical and, except for the final round
    if (l != numlevels-1) {
      res2 = exchange_shares(res1);
    }
  }

  // Return one share of the final result.
  // We need to call exchange again before using it in a subsequent operation.
  return res1;
}

// Asynchronous boolean equality (with eager rnum generation)
BShare eq_b_async(BShare x1, BShare x2, BShare y1, BShare y2) {

  // compute bitwise x^y^1
  BShare res1 = x1 ^ y1 ^ (~(BShare)0); // local share
  BShare res2 = x2 ^ y2 ^ (~(BShare)0); // remote share

  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);

  // The result is stored in the (numbits/2) rightmost bits of res1, res2
  for (int l=0; l<numlevels; l++) {
    res1 = eq_b_level2(numbits >> l, res1, res2);

    // Exchange results of logical and, except for the final round
    if (l != numlevels-1) {
      res2 = exchange_shares_async(res1);
    }
  }

  // Return one share of the final result.
  // We need to call exchange again before using it in a subsequent operation.
  return res1;
}

// Array-based boolean equality
void eq_b_array(BShare *x1, BShare *x2, BShare *y1, BShare *y2, long len,
                BShare *res) {
    BShare* res2 = (BShare *) malloc(len*sizeof(BShare)); // remote shares
    int numbits = sizeof(BShare) * 8;
    int rounds = (int) log2(numbits); // Number of rounds for oblivious equality
    BShare shifted_lbs, shifted_rbs;

    // compute bitwise x^y^1
    for (int i=0; i<len; i++) {
        res[i] = x1[i] ^ y1[i] ^ ~(BShare) 0; // local share;
        res2[i] = x2[i] ^ y2[i] ^ ~(BShare) 0; // remote share
    }

    for (int l=0; l<rounds-1; l++) {
        for (int i=0; i<len; i++) {
            shift_greater(res[i], res2[i], l, numbits, &shifted_lbs, &shifted_rbs);
            res[i] = and_b(res[i], res2[i], shifted_lbs, shifted_rbs, get_next_rb());
        }
        // Exchange shares
        exchange_shares_array(res, res2, len);
    }
    // Do final round without exchanging shares
    for (int i=0; i<len; i++) {
        shift_greater(res[i], res2[i], rounds-1, numbits, &shifted_lbs, &shifted_rbs);
        res[i] = and_b(res[i], res2[i], shifted_lbs, shifted_rbs, get_next_rb())
                 & (BShare) 1;
    }
    // The local share of the final result is stored in res.
    // Note: We need to call exchange again before using it in a subsequent operation.
    free(res2);
}

// array-based boolean equality with interleave between computation and communication
// per element
void eq_b_array_inter(BShare *x1, BShare *x2, BShare *y1, BShare *y2, long len,
                BShare *res) {

  BShare *res2 = (BShare *) malloc(len*sizeof(BShare)); // remote shares
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  MPI_Request r1, r2;
  int exchanges = 10;
  long batch_size = len/exchanges;

  // compute bitwise x^y^1
  for (long i=0; i<len; i++) {
    res[i] = x1[i] ^ y1[i] ^ (~(BShare)0); // local share;
    res2[i] = x2[i] ^ y2[i] ^ (~(BShare)0); // remote share
  }

  // The result is stored in the (numbits/2) rightmost bits of res, res2 elements
  for (int l=0; l<numlevels-1; l++) {

    // exchange 10 times per level
    for (long i=0; i<len; i++) {

      res[i] = eq_b_level2(numbits >> l, res[i], res2[i]);

      // first 9 exchanges
      if ( ((i+1)%batch_size)==0 ) {
        MPI_Irecv(&res2[i-(batch_size-1)], batch_size, MPI_LONG_LONG, get_succ(), XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r2);
        MPI_Isend(&res[i-(batch_size-1)], batch_size, MPI_LONG_LONG, get_pred(), XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
      }

      // last exchange
      if (i == len-1) {
        MPI_Wait(&r1, MPI_STATUS_IGNORE);
        MPI_Wait(&r2, MPI_STATUS_IGNORE);
        MPI_Irecv(&res2[len - batch_size], batch_size, MPI_LONG_LONG, get_succ(), XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r2);
        MPI_Isend(&res[len - batch_size], batch_size, MPI_LONG_LONG, get_pred(), XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
      }
    }
    // wait for last exchange before moving on the next level
    MPI_Wait(&r1, MPI_STATUS_IGNORE);
    MPI_Wait(&r2, MPI_STATUS_IGNORE);
  }

  // Last level (no exchange)
  for (long i=0; i<len; i++) {
    res[i] = eq_b_level2(numbits >> (numlevels-1), res[i], res2[i]);
  }

  // The local share of the final result is stored in res.
  // We need to call exchange again before using it in a subsequent operation.
  free(res2);
}


// array-based boolean equality with asynchronous exchange
// and interleaving between batches of elements
void eq_b_array_inter_batch(BShare *x1, BShare *x2, BShare *y1, BShare *y2, long len,
                BShare *res) {

  BShare *res2 = (BShare *) malloc(len*sizeof(BShare)); // remote shares
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  MPI_Request *r1 = (MPI_Request *) malloc(len*sizeof(MPI_Request));
  MPI_Request *r2 = (MPI_Request *) malloc(len*sizeof(MPI_Request));

  // compute bitwise x^y^1
  for (long i=0; i<len; i++) {
    res[i] = x1[i] ^ y1[i] ^ (~(BShare)0); // local share;
    res2[i] = x2[i] ^ y2[i] ^ (~(BShare)0); // remote share
  }

  for (int l=0; l<numlevels; l++) {

    // wait for results of previous level
    if (l > 0) {
      MPI_Waitall(len, r2, MPI_STATUSES_IGNORE);
      MPI_Waitall(len, r1, MPI_STATUSES_IGNORE);
    }
    for (long i=0; i<len; i++) {
      res[i] = eq_b_level2(numbits >> l, res[i], res2[i]);
      // exchange result for element i, level l
      // except for the final round
      if (l != numlevels-1) {
        MPI_Irecv(&res2[i], 1, MPI_LONG_LONG, get_succ(), XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r2[i]);
        MPI_Isend(&res[i], 1, MPI_LONG_LONG, get_pred(), XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1[i]);
      }
    }
  }
  // The local share of the final result is stored in res.
  // We need to call exchange again before using it in a subsequent operation.
  free(res2);
}

// same as eq_b_level but without rnums argument
BShare eq_b_level2(int numbits, BShare z1, BShare z2) {

  const BShare mask = 1;
  BShare res = 0;

  for (int i=0, j=0; i<numbits; i+=2, j++) {
    BShare bx1 = (z1 >> i) & mask; // bit at position i
    BShare by1 = (z1 >> (i+1)) & mask; // bit at position i+1
    BShare bx2 = (z2 >> i) & mask; // bit at position i of 2nd share
    BShare by2 = (z2 >> (i+1)) & mask; // bit at position i+1 of 2nd share

    // store the result (and's LSB) in the result's jth bit
    BShare out = and_b(bx1, bx2, by1, by2, get_next_rb());
    res |= (out & mask) << j;
  }
  return res;
}

// Computes x >? y using boolean shares in 1+logN communication rounds
// Avoids redundant computations by reusing results from previous rounds
BitShare greater(const BShare x1, const BShare x2,
                 const BShare y1, const BShare y2) {
    // Buffers
    BShare local[2];
    BShare remote[2];
    int len = sizeof(BShare)*8;   // Size of BShare in bits
    int rounds = (int) log2(len); // Number of rounds for oblivious inequality
    // Step 1: Compute l_i = x_i XOR y_i XOR 1, 1 <= i < BSHARE_LENGTH
    BShare local_x_xor_y = x1 ^ y1;
    BShare remote_x_xor_y = x2 ^ y2;
    BShare lbs = (local_x_xor_y ^ ~((BShare) 0));
    BShare rbs = (remote_x_xor_y ^ ~((BShare) 0));
    BShare shifted_lbs, shifted_rbs;
    // Compute first round
    shift_greater(lbs, rbs, 1, len, &shifted_lbs, &shifted_rbs);
    // Local bit shares of the 1st round
    local[0] = and_b(lbs, rbs,
                     shifted_lbs, shifted_rbs,
                     get_next_rb());
    // Set LSB to 'y_0 XOR 1' and compute the diagonal
    BShare l_diag = unset_lsbs(local_x_xor_y, 1);
    BShare r_diag = unset_lsbs(remote_x_xor_y, 1);
    l_diag |= (get_bit(y1, 0) ^ (BShare) 1);
    r_diag |= (get_bit(y2, 0) ^ (BShare) 1);
    // Local bit shares of the diagonal
    local[1] = and_b(l_diag, r_diag,
                     x1, x2,
                     get_next_rb());
    // Fetch remote bit shares
    exchange_shares_array(local, remote, 2);
    lbs = local[0]; rbs = remote[0];
    l_diag = local[1]; r_diag = remote[1];
    // Compute the remaining rounds
    for (int r=2; r<=rounds; r++) {
        shift_greater(lbs, rbs, r, len, &shifted_lbs, &shifted_rbs);
        lbs = and_b(lbs, rbs, shifted_lbs, shifted_rbs, get_next_rb());
        rbs = exchange_shares(lbs);
    }
    // Do a final AND with the diagonal
    lbs = and_b(set_bit(lbs>>1, len-1), set_bit(rbs>>1, len-1), // Shift to right and set MSB
                l_diag, r_diag,
                get_next_rb());
    rbs = exchange_shares(lbs);
    // XOR all bits and return local bit share
    BitShare res=0;
    for (int i=0; i<len; i++)
        res ^= get_bit(lbs, i);
    return res;
}

void greater_batch(const BShare *x1, const BShare *x2,
                  const BShare *y1, const BShare *y2,
                  int numElements, BitShare *res) {
    // Local and remote bits - We need two BShares per comparison
    unsigned long long *local_bits = (unsigned long long *) malloc(2 * numElements * sizeof(BShare));
    assert(local_bits!=NULL);
    unsigned long long *remote_bits = (unsigned long long *) malloc(2 * numElements * sizeof(BShare));
    assert(remote_bits!=NULL);
    int share_length = sizeof(BShare)*8; // The length of BShare in number of bits
    int rounds = (int) log2(share_length); // Number of rounds for oblivious inequality
    // Step 1: Compute l_i = x_i XOR y_i XOR 1, 1 <= i < BSHARE_LENGTH
    BShare local_x_xor_y, remote_x_xor_y, lbs, rbs, shifted_lbs, shifted_rbs, l_diag, r_diag;
    for (int i=0; i<numElements; i++) {
        local_x_xor_y = x1[i] ^ y1[i];
        remote_x_xor_y = x2[i] ^ y2[i];
        lbs = (local_x_xor_y ^ ~((BShare) 0));
        rbs = (remote_x_xor_y ^ ~((BShare) 0));
        // Compute first round
        shift_greater(lbs, rbs, 1, share_length, &shifted_lbs, &shifted_rbs);
        // Local bit shares of the 1st round
        local_bits[i] = and_b(lbs, rbs,
                              shifted_lbs, shifted_rbs,
                              get_next_rb());
        // Set LSB to 'y_0 XOR 1' and compute the diagonal
        l_diag = unset_lsbs(local_x_xor_y, 1);
        r_diag = unset_lsbs(remote_x_xor_y, 1);
        l_diag |= (get_bit(y1[i], 0) ^ (BShare) 1);
        r_diag |= (get_bit(y2[i], 0) ^ (BShare) 1);
        // Local bit shares of the diagonal (store it in the second half of the array)
        local_bits[numElements+i] = and_b(l_diag, r_diag,
                                          x1[i], x2[i],
                                          get_next_rb());
    }
    // Fetch remote bit shares
    exchange_shares_array_u(local_bits, remote_bits, numElements*2);
    // Compute the remaining rounds
    for (int r=2; r<=rounds; r++) {
        for (int i=0; i<numElements; i++) {
            lbs = local_bits[i];
            rbs = remote_bits[i];
            shift_greater(lbs, rbs, r, share_length, &shifted_lbs, &shifted_rbs);
            local_bits[i] = and_b(lbs, rbs, shifted_lbs, shifted_rbs, get_next_rb());
        }
        // Fetch remote shares
        exchange_shares_array_u(local_bits, remote_bits, numElements);
    }
    // Do a final AND with the diagonal
    for (int i=0; i<numElements; i++) {
        lbs = local_bits[i];
        rbs = remote_bits[i];
        l_diag = local_bits[numElements+i];
        r_diag = remote_bits[numElements+i];
        local_bits[i] = and_b(set_bit(lbs>>1, share_length-1), set_bit(rbs>>1, share_length-1),
                              l_diag, r_diag,
                              get_next_rb());
    }
    // Fetch remote shares
    exchange_shares_array_u(local_bits, remote_bits, numElements);
    // XOR all bits of each result
    for (int i=0; i<numElements; i++) {
        res[i]=0;
        for (int j=0; j<share_length; j++)
            res[i] ^= get_bit(local_bits[i], j);
    }
    // Free memory
    free(local_bits); free(remote_bits);
}

// Used by group_by_min_max_sel_odd_even()
void greater_batch2(BShare **c, int idx, int start, int numElements, 
                    int dist, BitShare *res) {
    // Local and remote bits - We need two BShares per comparison
    unsigned long long *local_bits = (unsigned long long *) malloc(2 * numElements * sizeof(BShare));
    assert(local_bits!=NULL);
    unsigned long long *remote_bits = (unsigned long long *) malloc(2 * numElements * sizeof(BShare));
    assert(remote_bits!=NULL);
    int share_length = sizeof(BShare)*8; // The length of BShare in number of bits
    int rounds = (int) log2(share_length); // Number of rounds for oblivious inequality
    // Step 1: Compute l_i = x_i XOR y_i XOR 1, 1 <= i < BSHARE_LENGTH
    BShare local_x_xor_y, remote_x_xor_y, lbs, rbs, shifted_lbs, shifted_rbs, l_diag, r_diag;
    for (int i=0, k=start; i<numElements; i++, k++) {
        local_x_xor_y = c[k][idx] ^ c[k+dist][idx];
        remote_x_xor_y = c[k][idx+1] ^ c[k+dist][idx+1];
        lbs = (local_x_xor_y ^ ~((BShare) 0));
        rbs = (remote_x_xor_y ^ ~((BShare) 0));
        // Compute first round
        shift_greater(lbs, rbs, 1, share_length, &shifted_lbs, &shifted_rbs);
        // Local bit shares of the 1st round
        local_bits[i] = and_b(lbs, rbs,
                              shifted_lbs, shifted_rbs,
                              get_next_rb());
        // Set LSB to 'y_0 XOR 1' and compute the diagonal
        l_diag = unset_lsbs(local_x_xor_y, 1);
        r_diag = unset_lsbs(remote_x_xor_y, 1);
        l_diag |= (get_bit(c[k+dist][idx], 0) ^ (BShare) 1);
        r_diag |= (get_bit(c[k+dist][idx+1], 0) ^ (BShare) 1);
        // Local bit shares of the diagonal (store it in the second half of the array)
        local_bits[numElements+i] = and_b(l_diag, r_diag,
                                          c[k][idx], c[k][idx+1],
                                          get_next_rb());
    }
    // Fetch remote bit shares
    exchange_shares_array_u(local_bits, remote_bits, numElements*2);
    // Compute the remaining rounds
    for (int r=2; r<=rounds; r++) {
        for (int i=0; i<numElements; i++) {
            lbs = local_bits[i];
            rbs = remote_bits[i];
            shift_greater(lbs, rbs, r, share_length, &shifted_lbs, &shifted_rbs);
            local_bits[i] = and_b(lbs, rbs, shifted_lbs, shifted_rbs, get_next_rb());
        }
        // Fetch remote shares
        exchange_shares_array_u(local_bits, remote_bits, numElements);
    }
    // Do a final AND with the diagonal
    for (int i=0; i<numElements; i++) {
        lbs = local_bits[i];
        rbs = remote_bits[i];
        l_diag = local_bits[numElements+i];
        r_diag = remote_bits[numElements+i];
        local_bits[i] = and_b(set_bit(lbs>>1, share_length-1), set_bit(rbs>>1, share_length-1),
                              l_diag, r_diag,
                              get_next_rb());
    }
    // Fetch remote shares
    exchange_shares_array_u(local_bits, remote_bits, numElements);
    // XOR all bits of each result
    for (int i=0; i<numElements; i++) {
        res[i]=0;
        for (int j=0; j<share_length; j++)
            res[i] ^= get_bit(local_bits[i], j);
    }
    // Free memory
    free(local_bits); free(remote_bits);
}

void geq_batch(const BShare *x1, const BShare *x2,
                 const BShare *y1, const BShare *y2,
                 int numElements, BitShare *res) {
    // Use greater_batch()
    greater_batch(y1, y2, x1, x2, numElements, res);
    // Negate result bits
    for (int i=0; i<numElements; i++)
        res[i] ^= (BitShare) 1;
}

// Computes x >=? y using boolean shares in 1+logN communication rounds
// Avoids redundant computations by reusing results from previous rounds
BitShare geq(const BShare x1, const BShare x2,
             const BShare y1, const BShare y2) {
    // x >=? y <==> NOT (x <? y)
    return greater(y1, y2, x1, x2) ^ (BitShare) 1;
}

// Compares and swaps two numbers x, y using their boolean shares
// Relies on greater() for number comparison, hence, it only works
// for numbers of the same sign
void cmp_swap(BShare *x1, BShare *x2, BShare *y1, BShare *y2,
              const BShare* rnums) {
  // Compute x > y
  BitShare b = greater(*x1, *x2, *y1, *y2);
  BShare bs1 = to_bshare(b);
  BShare bs2 = exchange_shares(bs1);
  BShare b1 = -bs1; // Set all bits equal to LSB of bs1
  BShare b2 = -bs2; // Set all bits equal to LSB of bs2
  BShare local[2];
  // Compute min = b * y + (1-b) * x
  local[0] = and_b(b1, b2, *y1, *y2, rnums[0]);
  local[0] ^= and_b(~b1, ~b2, *x1, *x2, rnums[1]);
  // Compute max = b * x + (1-b) * y
  local[1] = and_b(b1, b2, *x1, *x2, rnums[2]);
  local[1] ^= and_b(~b1, ~b2, *y1, *y2, rnums[3]);
  // Get remote shares from the other party
  BShare remote[2];
  exchange_shares_array(local, remote, 2);
  // Swap
  *x1 = local[0];
  *x2 = remote[0];
  *y1 = local[1];
  *y2 = remote[1];
}

// Generalized compare and swap
void cmp_swap_g(BShare* r1, BShare* r2, unsigned att_index1,
                unsigned att_index2, int num_elements, int asc) {
  // Compute x > y
  BitShare b = greater(r1[att_index1], r1[att_index1+1],
                       r2[att_index2], r2[att_index2+1]);
  BShare bs1 = to_bshare(b);
  BShare bs2 = exchange_shares(bs1);
  BShare b1 = -bs1; // Set all bits equal to LSB of bs1
  BShare b2 = -bs2; // Set all bits equal to LSB of bs2
  // Compute min, max for each pair of elements in the given arrays
  BShare local[num_elements];
  BShare r[2*num_elements];
  get_next_rb_array(r, 2*num_elements);
  for (int i=0, j=0; i<num_elements-1; i+=2, j+=4) {
    // Compute min = b * y + (1-b) * x
    local[i] = and_b(b1, b2, r2[i], r2[i+1], r[j]);
    local[i] ^= and_b(~b1, ~b2, r1[i], r1[i+1], r[j+1]);
    // Compute max = b * x + (1-b) * y
    local[i+1] = and_b(b1, b2, r1[i], r1[i+1], r[j+2]);
    local[i+1] ^= and_b(~b1, ~b2, r2[i], r2[i+1], r[j+3]);
  }
  // Get remote shares from the other party
  BShare remote[num_elements];
  exchange_shares_array(local, remote, num_elements);
  // Swap arrays
  int desc = !asc;
  for (int i=0; i<num_elements-1; i+=2) {
    r1[i] = local[i+desc];
    r1[i+1] = remote[i+desc];
    r2[i] = local[i+1-desc];
    r2[i+1] = remote[i+1-desc];
  }
}

// Adds to numbers using their boolean shares
// This function represents a single bit using
// a boolean share. As a result, it requires exchanging two BShares for each
// bitwise logical AND, whereas we should only exchange two bits
void boolean_addition(BShare x1, BShare x2, BShare y1, BShare y2,
                      BShare *z1, BShare *z2,
                      BShare *rnums) {

    BShare c1[3], c2[3], xor_b=0, carry1=0, carry2=0, mask=1;
    // Add bits ingoring carries
    *z1 = x1 ^ y1;
    *z2 = x2 ^ y2;
    // Compute carry at position 0 (x_0 AND y_0)
    carry1 = and_b(get_bit(x1, 0), get_bit(x2, 0),
                   get_bit(y1, 0), get_bit(y2, 0), rnums[0]) & mask;
    carry2 = exchange_shares(carry1) & mask;
    // Compute carry at position i, ignoring the possible overflow carry
    for (int i=1; i<(sizeof(BShare)*8) - 1; i++) {
      // x_i AND y_i
      c1[0] = and_b(get_bit(x1, i), get_bit(x2, i),
                    get_bit(y1, i), get_bit(y2, i),
                    rnums[i]);
      // x_i AND c_{i-1}
      c1[1] = and_b(get_bit(x1, i), get_bit(x2, i),
                    get_bit(carry1, i-1),
                    get_bit(carry2, i-1),
                    rnums[i+1]);
      // y_i AND c_{i-1}
      c1[2] =  and_b(get_bit(y1, i), get_bit(y2, i),
                     get_bit(carry1, i-1),
                     get_bit(carry2, i-1),
                     rnums[i+2]);
      // Get carry shares from the other party
      exchange_shares_array(c1, c2, 3);
      // Store locally computed bit in carry1
      xor_b = (c1[0] ^ c1[1] ^ c1[2]);
      carry1 |= ((xor_b & mask) << i);
      // Store other party's bit share in carry2
      xor_b = (c2[0] ^ c2[1] ^c2[2]);
      carry2 |= ((xor_b & mask) << i);
    }
    // Add carries
    *z1 ^= (carry1 << 1);
    *z2 ^= (carry2 << 1);
}

void boolean_addition_batch(BShare *x1, BShare *x2, BShare *y1, BShare *y2,
                       BShare *res, int numElements) {

  BShare xor_b=0;
  BShare **c1 = allocate_2D_table(numElements, 3);
  BShare **c2 = allocate_2D_table(numElements, 3);

  BShare *carry1 = (BShare *) malloc(numElements * sizeof(BShare));
  BShare *carry2 = (BShare *) malloc(numElements * sizeof(BShare));
  BShare mask = 1;

  // Add bits ingoring carries
  for (int i=0; i<numElements; i++) {
    res[i] = x1[i]^y1[i];
    // Compute carry at position 0 (x_0 AND y_0)
    carry1[i] = and_b(get_bit(x1[i], 0), get_bit(x2[i], 0),
                   get_bit(y1[i], 0), get_bit(y2[i], 0),
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
      c1[i][0] = and_b(get_bit(x1[i], j), get_bit(x2[i], j),
                    get_bit(y1[i], j), get_bit(y2[i], j),
                    get_next_rb());
      // x_j AND c_{j-1}
      c1[i][1] = and_b(get_bit(x1[i], j), get_bit(x2[i], j),
                    get_bit(carry1[i], j-1),
                    get_bit(carry2[i], j-1),
                    get_next_rb());
      // y_i AND c_{i-1}
      c1[i][2] =  and_b(get_bit(y1[i], j), get_bit(y2[i], j),
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

void boolean_addition_batch2(BShare **c, BShare *res1, 
                             int start, int dist, int numElements,
                             int idx) {
  
  BShare xor_b=0;
  BShare **c1 = allocate_2D_table(numElements, 3);
  BShare **c2 = allocate_2D_table(numElements, 3);

  BShare *carry1 = (BShare *) malloc(numElements * sizeof(BShare));
  BShare *carry2 = (BShare *) malloc(numElements * sizeof(BShare));
  BShare mask = 1;

  // Add bits ingoring carries
  for (int i=start, k=0; k<numElements; i++, k++) {
    res1[k] = c[i][idx]^c[i+dist][idx];
    // Compute carry at position 0 (x_0 AND y_0)
    carry1[k] = and_b(get_bit(c[i][idx], 0), get_bit(c[i][idx+1], 0),
                   get_bit(c[i+dist][idx], 0), get_bit(c[i+dist][idx+1], 0),
                   get_next_rb()) & mask;
  }

  exchange_shares_array(carry1, carry2, numElements);
  for (int i=0; i<numElements; i++) {
    carry2[i] &= mask;
  }

  int bsize = sizeof(BShare)*8;

  // Compute carry at position j, ignoring the possible overflow carry
  for (int j=1; j<bsize - 1; j++) {
    for (int i=start, k=0; k<numElements; i++, k++) {
      // x_j AND y_j
      c1[k][0] = and_b(get_bit(c[i][idx], j), get_bit(c[i][idx+1], j),
                    get_bit(c[i+dist][idx], j), get_bit(c[i+dist][idx+1], j),
                    get_next_rb());
      // x_j AND c_{j-1}
      c1[k][1] = and_b(get_bit(c[i][idx], j), get_bit(c[i][idx+1], j),
                    get_bit(carry1[k], j-1),
                    get_bit(carry2[k], j-1),
                    get_next_rb());
      // y_i AND c_{i-1}
      c1[k][2] =  and_b(get_bit(c[i+dist][idx], j), get_bit(c[i+dist][idx+1], j),
                     get_bit(carry1[k], j-1),
                     get_bit(carry2[k], j-1),
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
    res1[i] ^= (carry1[i] << 1);
  }
  free(c1); free(c2); free(carry1); free(carry2);
}

// Less than 0
BShare ltz_b(BShare x) {
  return (unsigned long long) x >> (sizeof(x)*8 -1);
}

// Less than 0 for an array of elements
void ltz_b_array(const BShare *x, int len, BShare *res) {
  for (int i=0; i<len; i++) {
    res[i] = ltz_b(x[i]);
  }
}

// Conversion of a single-bit BShare to an arithmetic share
AShare convert_single_bit(BShare bit, AShare ra, BShare rb) {
  // reveal bit^rb
  Data z = reveal_b((bit&1)^(rb&1));
  // return [ra] * (1-2*z) + z
  // Only rank 2 needs to subtract 1 if z=1, so we multiply by rank mod 2
  ra = ra * (1-2*z) + z*(get_rank()%2);
  return ra;
}

// Conversion of a single-bit BShare array to an array of arithmetic shares
void convert_single_bit_array(BShare *bit, AShare *ra, BShare *rb, int len,
                               AShare *res) {
  int i;
  BShare* z = (BShare*) malloc(len*sizeof(BShare));
  assert(z!=NULL);
  for (i=0; i<len; i++) {
    z[i] = (bit[i]&1)^(rb[i]&1);
  }
  reveal_b_array_async(z, len);
  // return [ra] * (1-2*z) + z
  // Only rank 2 needs to subtract 1 if z=1, so we multiply by rank mod 2
  for (i=0; i<len; i++) {
    res[i] = ra[i] * (1-2*z[i]) + z[i]*(get_rank()%2);
  }
  free(z);
}

// Used by cmp_swap_batch() to get the next valid index
inline PRIVATE int get_next_index(int pos, int area, int comp_per_box) {
  int box_start = (pos / area) * area;
  if ((pos + 1) >= (box_start + comp_per_box)) {
    return box_start + area;
  }
  return pos + 1;
}

// converts an arithmetic share to binary
void convert_a_to_b_array(AShare *xa1, AShare *xa2, BShare *xb1, BShare *xb2, int len) {

  MPI_Request r1, r2;
  BShare *w1 = (BShare *) malloc(len*sizeof(BShare)); // local share of x3
  assert(w1!=NULL);
  BShare *w2 = (BShare *) malloc(len*sizeof(BShare)); // remote share of x3
  assert(w2!=NULL);
  BShare *r_temp = (BShare *) malloc(len*sizeof(BShare)); // random share
  assert(r_temp!=NULL);
  BShare *r_temp2 = (BShare *) malloc(len*sizeof(BShare)); // random share
  assert(r_temp2!=NULL);

  if (get_rank() == 0) {
    // generate bool shares of x1+x2 (xa1 + xa2)
    BShare *z13 = (BShare *) malloc(len*sizeof(BShare));
    assert(z13!=NULL);
    BShare *z12 = (BShare *) malloc(len*sizeof(BShare));
    assert(z12!=NULL);
    for (int i=0; i<len; i++) {
      xa2[i] += xa1[i]; // xa1 + xa2
    }
    for (int i=0; i<len; i++) {
      generate_bool_share(xa2[i], &xa1[i], &z12[i], &z13[i]);
      // local share of xa1+xa2 now in xa1
    }
    // distribute shares to P2, P3
    MPI_Isend(z12, len, MPI_LONG_LONG, 1, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
    MPI_Isend(z13, len, MPI_LONG_LONG, 2, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r2);
    MPI_Wait(&r1, MPI_STATUS_IGNORE);
    MPI_Wait(&r2, MPI_STATUS_IGNORE);
    free(z12); free(z13);

    // receive share of x3 from P3
    MPI_Irecv(w1, len, MPI_LONG_LONG, 2, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
    MPI_Wait(&r1, MPI_STATUS_IGNORE);

    // generate pairs of random binary shares (R1)
    get_next_rb_pair_array(xb2, xb1, len);
    // receive share from P2 and compute xb1
    MPI_Irecv(r_temp, len, MPI_LONG_LONG, 1, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
    MPI_Wait(&r1, MPI_STATUS_IGNORE);
    // compute xb1
    for (int i=0; i<len; i++) {
      xb1[i] ^= r_temp[i];
      xb1[i] ^= xb2[i]; // R1
    }
    // xb2 contains the local share of R1
    // generate pairs of random binary shares (R2)
    get_next_rb_pair_array(r_temp, r_temp2, len);
    // r_temp contains the local share of R2
  }
  else if (get_rank() == 1) { //P2
    // receive share of x1+x2 frm P1
    MPI_Irecv(xa1, len, MPI_LONG_LONG, 0, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
    // receive share of x3
    MPI_Irecv(w1, len, MPI_LONG_LONG, 2, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r2);
    MPI_Wait(&r1, MPI_STATUS_IGNORE);
    MPI_Wait(&r2, MPI_STATUS_IGNORE);

    // generate pairs of random binary shares (R1)
    get_next_rb_pair_array(xb2, xb1, len);
    // send local to P1
    MPI_Isend(xb2, len, MPI_LONG_LONG, 0, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
    MPI_Wait(&r1, MPI_STATUS_IGNORE);
    // xb2 contains the local share of R1
    // generate pairs of random binary shares (R2)
    get_next_rb_pair_array(r_temp, xb1, len);
    // receive share from P3 and compute xb1
    MPI_Irecv(r_temp2, len, MPI_LONG_LONG, 2, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
    MPI_Wait(&r1, MPI_STATUS_IGNORE);
    // compute xb1
    for (int i=0; i<len; i++) {
      xb1[i] ^= r_temp2[i];
      xb1[i] ^= r_temp[i]; // R2
    }
    // r_temp contains the local share of R2
  }
  else { //P3
    // generate bool shares of x3 (xa1)
    BShare *w13 = (BShare *) malloc(len*sizeof(BShare));
    assert(w13!=NULL);

    for (int i=0; i<len; i++) {
      generate_bool_share(xa1[i], &w13[i], &w2[i], &w1[i]);
    }

    // receive share of x1+x2 frm P1
    MPI_Irecv(xa1, len, MPI_LONG_LONG, 0, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
    MPI_Wait(&r1, MPI_STATUS_IGNORE);

    // distribute shares to P1, P2
    MPI_Isend(w13, len, MPI_LONG_LONG, 0, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
    MPI_Isend(w2, len, MPI_LONG_LONG, 1, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r2);

    MPI_Wait(&r1, MPI_STATUS_IGNORE);
    MPI_Wait(&r2, MPI_STATUS_IGNORE);
    free(w13);

    // generate pairs of random binary shares (R1)
    get_next_rb_pair_array(xb2, xb1, len);
    // xb2 contains the local share of R1
    // generate pairs of random binary shares (R2)
    get_next_rb_pair_array(r_temp, xb1, len);
    // send local to P2
    MPI_Isend(r_temp, len, MPI_LONG_LONG, 1, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
    MPI_Wait(&r1, MPI_STATUS_IGNORE);
    // r_temp contains the local share of R2
  }

  /***** all parties *****/
  // get remote shares z2, w2
  exchange_shares_array(xa1, xa2, len);
  exchange_shares_array(w1, w2, len);
  // ripple-carry-adder
  boolean_addition_batch(xa1, xa2, w1, w2, r_temp2, len); // r_temp2: result of RCA
  free(w2);

  // share of third share
  for (int i=0; i<len; i++) {
    w1[i] = r_temp2[i] ^ xb2[i] ^ r_temp[i];
  }
  free(r_temp); free(r_temp2);

  // reveal y to P3
  if (get_rank() == 0) {
      MPI_Isend(w1, len, MPI_LONG_LONG, 2, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
      MPI_Wait(&r1, MPI_STATUS_IGNORE);
  } else if (get_rank() == 1) { // P2
      MPI_Isend(w1, len, MPI_LONG_LONG, 2, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
      MPI_Wait(&r1, MPI_STATUS_IGNORE);
  } else { // P3
      MPI_Irecv(xb1, len, MPI_LONG_LONG, 0, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r1);
      MPI_Irecv(xa1, len, MPI_LONG_LONG, 1, XCHANGE_MSG_TAG, MPI_COMM_WORLD, &r2);

      MPI_Wait(&r1, MPI_STATUS_IGNORE);
      MPI_Wait(&r2, MPI_STATUS_IGNORE);
      for (int i=0; i<len; i++) {
        xb1[i] ^= xa1[i];
        xb1[i] ^= w1[i];
      }
    }
  /***** all parties *****/
  // exchange xb1, xb2
  exchange_shares_array(xb1, xb2, len);
  free(w1);
}

// Used by bitonic_sort_batch()
// Applies each round of cmp_swap_g() for all 'num_comparisons'
int cmp_swap_batch(int first_index, BShare** rows, int length, int row_length,
                   unsigned* att_indices, int num_attributes, int phase,
                   int column, bool* asc, int num_comparisons) {

    // Number of 'boxes' per column of bitonic sort
    long long boxes = length / (((long long) 1) << (phase + 1 - column));
    // Number of comparisons per 'box'
    int comp_per_box = length / (2 * boxes);
    // The 'area' of the 'box'
    int area = 2 * comp_per_box;

    int numbits = sizeof(BShare) * 8;
    int numlevels = log2(numbits);

    // We keep 'num_attributes' bits per comparison
    BitShare** bs1 = allocate_2D_bit_table(num_comparisons, num_attributes);
    BitShare** bs2 = allocate_2D_bit_table(num_comparisons, num_attributes);

    // We keep 'num_attributes' bits per comparison
    BitShare** eq_bs1 = NULL;
    BitShare** eq_bs2 = NULL;
    if (num_attributes>1) { // Only if we have more than one sort attributes
        eq_bs1 = allocate_2D_bit_table(num_comparisons, num_attributes);
        eq_bs2 = allocate_2D_bit_table(num_comparisons, num_attributes);
    }

    // We perform 'num_comparisons' comparisons per 'column'
    unsigned long long *local = (unsigned long long *) malloc(2 * num_comparisons * sizeof(BShare));
    assert(local!=NULL);
    unsigned long long *remote = (unsigned long long *) malloc(2 * num_comparisons * sizeof(BShare));
    assert(remote!=NULL);

    int share_len = sizeof(BShare)*8;   // Size of BShare in bits
    int rounds = (int) log2(share_len); // Number of rounds for oblivious inequality

    int done, pos;
    // Apply greater_b2() and eq_b_array() for each sorting attribute
    for (int att=0; att<num_attributes; att++) {
        // The sorting attribute
        int att_index = att_indices[att];

        // Distance between two compared elements
        int d = 1 << (phase-column);

        // 1. Do 1st round of greater2() for the whole batch
        done = 0;           // Number of comparisons done
        pos = first_index;  // Index of the first element in comparison
        BShare local_x_xor_y, remote_x_xor_y,
                lbs, rbs,
                shifted_lbs, shifted_rbs,
                l_diag, r_diag;
        while ((pos < length) & (done != num_comparisons)) {
            local_x_xor_y = rows[pos][att_index] ^ rows[pos|d][att_index];
            remote_x_xor_y = rows[pos][att_index+1] ^ rows[pos|d][att_index+1];
            lbs = (local_x_xor_y ^ ~((BShare) 0));
            rbs = (remote_x_xor_y ^ ~((BShare) 0));
            // Compute first round
            shift_greater(lbs, rbs, 1, share_len, &shifted_lbs, &shifted_rbs);
            // Local bit shares of the 1st round
            local[done] = and_b(lbs, rbs,
                                shifted_lbs, shifted_rbs,
                                get_next_rb());
            // Set LSB to 'y_0 XOR 1' (or to 'x_0 XOR 1', depending on ASC/DESC) and compute the diagonal
            l_diag = unset_lsbs(local_x_xor_y, 1);
            r_diag = unset_lsbs(remote_x_xor_y, 1);
            l_diag |= (get_bit(rows[pos|(d*asc[att])][att_index], 0) ^ (BShare) 1);
            r_diag |= (get_bit(rows[pos|(d*asc[att])][att_index+1], 0) ^ (BShare) 1);
            // Local bit shares of the diagonal (store it in the second half of the array)
            local[num_comparisons+done] = and_b(l_diag, r_diag,
                                                rows[pos|(d*!asc[att])][att_index],
                                                rows[pos|(d*!asc[att])][att_index+1],
                                                get_next_rb());
            pos = get_next_index(pos, area, comp_per_box);
            done++;
        }
        // Fetch remote bits
        exchange_shares_array_u(local, remote, num_comparisons*2); // 1 round
        // Compute the remaining rounds
        for (int r=2; r<=rounds; r++) {
            for (int i=0; i<done; i++) {
                lbs = local[i];
                rbs = remote[i];
                shift_greater(lbs, rbs, r, share_len, &shifted_lbs, &shifted_rbs);
                local[i] = and_b(lbs, rbs, shifted_lbs, shifted_rbs, get_next_rb());
            }
            // Fetch remote shares
            exchange_shares_array_u(local, remote, done); // 1 round
        }
        // Do a final AND with the diagonal
        for (int i=0; i<done; i++) {
            lbs = local[i];
            rbs = remote[i];
            l_diag = local[num_comparisons+i];
            r_diag = remote[num_comparisons+i];
            local[i] = and_b(set_bit(lbs>>1, share_len-1), set_bit(rbs>>1, share_len-1),
                             l_diag, r_diag,
                             get_next_rb());
        }
        // Fetch remote shares
        exchange_shares_array_u(local, remote, done);
        // XOR all bits of each result
        for (int i=0; i<done; i++) {
            bs1[i][att]=0;
            for (int j=0; j<share_len; j++)
                bs1[i][att] ^= get_bit(local[i], j);
        }

        // Now compute equality bits (another logL rounds)
        if ( (num_attributes>1) && (att<num_attributes-1) ) {
            // Only if there are more than one sort attributes
            BShare* res = (BShare*) malloc(num_comparisons*sizeof(BShare));
            assert(res!=NULL);
            BShare* res2 = (BShare*) malloc(num_comparisons*sizeof(BShare));
            assert(res2!=NULL);
            done = 0;           // Number of equality comparisons done
            pos = first_index;  // Index of the first element in comparison
            while ((pos < length) & (done != num_comparisons)) {
                BShare x1 = rows[pos][att_index];
                BShare x2 = rows[pos][att_index+1];
                BShare y1 = rows[pos|d][att_index];
                BShare y2 = rows[pos|d][att_index+1];
                res[done] = x1 ^ y1 ^ (~(BShare)0); // local share;
                res2[done] = x2 ^ y2 ^ (~(BShare)0); // remote share

                pos = get_next_index(pos, area, comp_per_box);
                done++;
            }
            // The result is stored in the (numbits/2) rightmost bits of res, res2 elements
            for (int l=0; l<numlevels; l++) {
                done = 0;           // Number of equality comparisons done
                pos = first_index;  // Index of the first element in comparison
                while ((pos < length) & (done != num_comparisons)) {
                    res[done] = eq_b_level2(numbits >> l, res[done], res2[done]);

                    pos = get_next_index(pos, area, comp_per_box);
                    done++;
                }
                // Exchange results of logical and, except for the final round
                if (l != numlevels-1) {
                    exchange_shares_array(res, res2, done);
                }
            }
            BitShare mask = 1;
            // Copy equality bits for current attribute
            for (int c=0; c<done; c++) {
                eq_bs1[c][att] = (res[c] & mask);
            }
            free(res); free(res2);
        }
    }
    free(local); free(remote);

    // Get remote bit shares
    exchange_bit_shares_array(&bs1[0][0], &bs2[0][0],
                              num_comparisons*num_attributes);   // 1 round
    if (num_attributes>1){
        exchange_bit_shares_array(&eq_bs1[0][0], &eq_bs2[0][0],
                                  num_comparisons*num_attributes);   // 1 round
    }

    // Compute composite bit for each comparison
    if (num_attributes==2) {
        compute_composite_2(bs1, bs2, eq_bs1, eq_bs2, num_comparisons,
                            num_attributes);
    }
    else if (num_attributes==3){ // num_attributes==3
        compute_composite_3(bs1, bs2, eq_bs1, eq_bs2, num_comparisons,
                            num_attributes);
    }

    if (num_attributes>1) {
        free(eq_bs1); free(eq_bs2);
    }

    // 4. Multiplexing
    BShare **local_rows = allocate_2D_table(num_comparisons, row_length);
    BShare **remote_rows = allocate_2D_table(num_comparisons, row_length);
    // Distance between two compared elements
    int d = 1 << (phase-column);
    done = 0;
    pos = first_index;
    while ((pos < length) & (done != num_comparisons)) {
        BShare b1 = - (BShare) bs1[done][0]; // Set all bits equal to LSB
        BShare b2 = - (BShare) bs2[done][0]; // Set all bits equal to LSB
        // Compute min, max for each pair of elements in the given arrays
        BShare* r1 = rows[pos];
        BShare* r2 = rows[pos|d];
        BShare r[2*row_length];
        // We need 4*row_length/2 = 2*row_length random numbers
        get_next_rb_array(r, 2*row_length);
        // For each row element
        for (int j=0, k=0; j<row_length-1; j+=2, k+=4) {
            // Compute min = b * y + (1-b) * x
            local_rows[done][j] = and_b(b1, b2, r2[j], r2[j+1], r[k]);
            local_rows[done][j] ^= and_b(~b1, ~b2, r1[j], r1[j+1], r[k+1]);
            // Compute max = b * x + (1-b) * y
            local_rows[done][j+1] = and_b(b1, b2, r1[j], r1[j+1], r[k+2]);
            local_rows[done][j+1] ^= and_b(~b1, ~b2, r2[j], r2[j+1], r[k+3]);
        }
        pos = get_next_index(pos, area, comp_per_box);
        done++;
    }
    // Get remote shares from the other party  -- 1 round
    exchange_shares_array(&local_rows[0][0], &remote_rows[0][0],
                          num_comparisons*row_length);

    // 5. Update table rows
    done = 0;
    pos = first_index;
    while ((pos < length) & (done != num_comparisons)) {
        // up=true means that max should be placed at the second slot,
        // otherwise at the first one
        int up = ((pos >> phase) & 2) == 0;
        int i = up ^ 1;
        // For each element in the row
        for (int j=0; j<row_length-1; j+=2) {
            rows[pos][j] = local_rows[done][j+i];
            rows[pos][j+1] = remote_rows[done][j+i];
            rows[pos|d][j] = local_rows[done][j+1-i];
            rows[pos|d][j+1] = remote_rows[done][j+1-i];
        }
        pos = get_next_index(pos, area, comp_per_box);
        done++;
    }

    // Free memory
    free(local_rows); free(remote_rows);
    free(bs1); free(bs2);

    // Return the last index checked for next call
    return pos;
}

// Compute composite b = b^1_g OR (b^1_e AND b^2_g) OR
//                       (b^1_e AND b^2_e AND b^3_g) OR
//                       ...
//                       (b^1_e AND b^2_e AND ... AND b^(n-1)_e AND b^n_g)
// This can be done in O(logn) rounds, where n is the number of sort attributes

// Computes composite b = b^1_g OR (b^1_e AND b^2_g)
// Composite bit shares are stored in bg1[0], bg2[0]
// Requires 2 communication rounds in total (independent from 'num_rows')
PRIVATE void compute_composite_2(BitShare** bg1, BitShare** bg2,
                                 BitShare** be1, BitShare** be2,
                                 int num_rows, int num_cols) {
  assert(num_cols==2);
  BShare mask=1;
  // 1st round: for each row
  for (int i=0; i<num_rows; i++) {
    // Compute b^1_e AND b^2_g and store result in be
    be1[i][0] = and_b(be1[i][0], be2[i][0],
                      bg1[i][1], bg2[i][1], get_next_rb())
                    & mask;
  }
  // Get remote shares for all bits -- 1 round
  exchange_bit_shares_array(&be1[0][0], &be2[0][0], num_rows*num_cols);

  // 2nd round: for each row
  for (int i=0; i<num_rows; i++) {
    // Compute NOT ( NOT(b^1_g) AND NOT(b^1_e AND b^2_g) )
    bg1[i][0] = and_b(bg1[i][0] ^ mask, bg2[i][0] ^ mask,
                      be1[i][0] ^ mask, be2[i][0] ^ mask, get_next_rb())
                    & mask;
    bg1[i][0] ^= mask;
  }
  // Fetch remote shares -- 1 round
  exchange_bit_shares_array(&bg1[0][0], &bg2[0][0], num_rows*num_cols);
}

// Computes composite b = b^1_g OR (b^1_e AND b^2_g AND NOT(b^2_e)) OR
//                        (b^1_e AND b^2_e AND b^3_g)
// Composite bit shares are stored in bg1[0], bg2[0]
// Requires 4 communication rounds in total (independent from 'num_rows')
PRIVATE void compute_composite_3(BitShare** bg1, BitShare** bg2,
                                 BitShare** be1, BitShare** be2,
                                 int num_rows, int num_cols) {
  assert(num_cols==3);
  BShare mask=1;
  BitShare** b1 = allocate_2D_bit_table(num_rows, 2);
  BitShare** b2 = allocate_2D_bit_table(num_rows, 2);

  // 1st round: for each row
  for (int i=0; i<num_rows; i++) {
    // Compute b^1_e AND b^2_g
    b1[i][0] = and_b(be1[i][0], be2[i][0], bg1[i][1], bg2[i][1], get_next_rb())
                  & mask;
    // Compute b^1_e AND b^2_e
    b1[i][1]= and_b(be1[i][0], be2[i][0], be1[i][1], be2[i][1], get_next_rb())
                  & mask;
  }
  // Fetch remote shares for all rows
  exchange_bit_shares_array(&b1[0][0], &b2[0][0], num_rows*2);    // 1 round

  // 2nd round: for each row
  for (int i=0; i<num_rows; i++) {
    // Compute (b^1_e AND b^2_g) AND NOT(b^2_e)
    b1[i][0] = and_b(b1[i][0], b2[i][0],
                     be1[i][1] ^ mask, be2[i][1] ^ mask, get_next_rb())
                    & mask;
    // Compute (b^1_e AND b^2_e) AND b^3_g
    b1[i][1] = and_b(b1[i][1], b2[i][1], bg1[i][2], bg2[i][2], get_next_rb())
                    & mask;
  }
  // Fetch remote shares for all rows
  exchange_bit_shares_array(&b1[0][0], &b2[0][0], num_rows*2);    // 1 round

  // 3rd round: for each row
  for (int i=0; i<num_rows; i++) {
    // Compute ( NOT(b^1_g) AND NOT(b^1_e AND NOT(b^2_e) AND b^2_g) )
    b1[i][0] = and_b(bg1[i][0] ^ mask, bg2[i][0] ^ mask,
                     b1[i][0] ^ mask, b2[i][0] ^ mask, get_next_rb())
                    & mask;
  }
  // Fetch remote shares for all rows
  exchange_bit_shares_array(&b1[0][0], &b2[0][0], num_rows*2);    // 1 round

  // 4th round: for each row
  for (int i=0; i<num_rows; i++) {
    // Compute ( NOT(b^1_g) AND NOT(b^1_e AND NOT(b^2_e) AND b^2_g) )
    //           AND NOT(b^1_e AND b^2_e AND b^3_g)
    b1[i][0] = and_b(b1[i][0], b2[i][0],
                     b1[i][1] ^ mask, b2[i][1] ^ mask, get_next_rb())
                    & mask;
  }
  // Fetch remote shares for all rows
  exchange_bit_shares_array(&b1[0][0], &b2[0][0], num_rows*2);    // 1 round
  // Store back composite bit shares
  for (int i=0; i<num_rows; i++) {
    bg1[i][0] = b1[i][0] ^ mask;
    bg2[i][0] = b2[i][0] ^ mask;
  }
  free(b1); free(b2);
}
