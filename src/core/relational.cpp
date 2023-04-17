#include "../../include/core/relational.h"

#define PRIVATE static

// Private function declarations
PRIVATE void distinct_batch_incr(BShareTable*, int, int, unsigned, BShare*);
PRIVATE void in_level(int, BShare*, BShare*);
PRIVATE void join_eq_b(BShareTable input1, BShareTable input2,
                       int leftcol, int rightcol, BShare result[]);
PRIVATE void select_eq(AShareTable input, int leftcol, AShare c1, AShare c2,
                       AShare result[]);
PRIVATE void select_eq_b(BShareTable input, int leftcol, int rightcol,
                         BShare result[]);
PRIVATE void select_geq_b(BShareTable input, int leftcol, BShare result[]);
PRIVATE void select_greater_b(BShareTable input, int leftcol, BShare result[]);
PRIVATE void select_greater_batch(BShareTable input, int leftcol, int rightcol,
                                  BShare result[]);
PRIVATE void shift_greater(const BShare, const BShare, int, int, BShare*, BShare*);

inline PRIVATE void shift_greater(const BShare x1, const BShare x2, int r, int len, BShare *shifted_1, BShare *shifted_2) {
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

// Selection on arithmetic shares
void select_a(AShareTable input, Predicate p, AShare c1, AShare c2, AShare result[]) {
  switch (p.operation) {
      case EQ:
          //printf("Equality predicate\n");
          select_eq(input, p.leftcol, c1, c2, result);
          break;
      default:
          printf("Illegal operation. Only equality is supported for the moment.\n");
  }
}

// Internal equality select: right must be a pointer to a pair of constants
PRIVATE void select_eq(AShareTable input, int leftcol, AShare c1, AShare c2, AShare result[]) {
  for (int i = 0; i < input.numRows ;i++) {

    // generate w and r
    WSharePair w = get_next_w();
    AShare r = get_next_r();

    // compute the equality result
    result[i] = eq(input.content[i][leftcol], // 1st share of left att
                    input.content[i][leftcol + 1], // 2nd share of left att
                    c1, // 1st share of right att (constant)
                    c2, // 2nd share left att (constant)
                    w.first, // 1st share of random w
                    w.second, // 2nd share of random w
                    r); // 1st share of random r
  }
}

void and_b_table(BShareTable input, int leftcol, int rightcol, int size, BShare result[]){
  for (int i = 0; i < size; i++){
    result[i] = and_b(input.content[i][leftcol],
                      input.content[i][leftcol + 1],
                      input.content[i][rightcol],
                      input.content[i][rightcol + 1], get_next_rb());
  }
}

// Selection on boolean shares
void select_b(const BShareTable input, Predicate_B p, BShare result[]) {
  switch (p.operation) {
    case EQ:
      select_eq_b(input, p.leftcol, p.rightcol, result);
      break;
    case GT:
      select_greater_b(input, p.leftcol, result);
      break;
    case GEQ:
      select_geq_b(input, p.leftcol, result);
      break;
    case GC:
      select_greater_batch_const(input, p.leftcol, p.cs1, p.cs2, result);
      break;
    case GR:
      select_greater_batch(input, p.leftcol, p.rightcol, result);
      break;
    case EQC:
      select_eq_batch_const(input, p.leftcol, p.cs1, p.cs2, result);
      break;
    default:
          printf("Illegal operation. Predicate %u not supported.\n", p.operation);
  }
}

/**
 *  Internal inequality select for boolean shares, i.e. att > c, where c is a public constant.
 *  This is copmputed as (c-att) < 0 and col points to the difference shares of the
 *  attribute to be tested:
 **/
PRIVATE void select_greater_b(BShareTable input, int leftcol, BShare result[]) {
  for (int i = 0; i < input.numRows; i++) {
    result[i] = ltz_b(input.content[i][leftcol]);
  }
}

/**
 *  Internal inequality select for boolean shares, i.e. att >= c, where c is a public constant.
 *  This is copmputed as ~((att - c) < 0) and col points to the difference shares of the
 *  attribute to be tested:
 **/
PRIVATE void select_geq_b(BShareTable input, int leftcol, BShare result[]) {
  for (int i = 0; i < input.numRows; i++) {
    result[i] = (ltz_b(input.content[i][leftcol]) & 1) ^ 1;
  }
}

/**
 *  Internal equality select for boolean shares.
 *  leftcol and rightcol point to the difference shares of the attribute to be tested:
 *  leftcol: att-c, rightcol: c-att
 **/
PRIVATE void select_eq_b(BShareTable input, int leftcol, int rightcol, BShare result[]) {
  // r = (c-att < 0) ^ (att-c < 0) ^ 1
  BShare z1, z2;
  for (int i = 0; i < input.numRows; i++) {
    z1 = ltz_b(input.content[i][leftcol]);
    z2 = ltz_b(input.content[i][rightcol]);
    result[i] = z1^z2^1;
  }
}

/**
 *  Internal inequality select for boolean shares.
 *  leftcol and rightcol point to two table columns
 **/
PRIVATE void select_greater_batch(BShareTable input, int leftcol, int rightcol,
                                  BShare result[]) {
    int numElements = input.numRows;
    BShare** c = input.content;
    int share_length=sizeof(BShare)*8; // The length of BShare in number of bits
    // Local and remote bits - We need two BShares per comparison
    unsigned long long *local_bits = (unsigned long long *) malloc(2 * numElements * sizeof(BShare));
    assert(local_bits!=NULL);
    unsigned long long *remote_bits = (unsigned long long *) malloc(2 * numElements * sizeof(BShare));
    assert(remote_bits!=NULL);
    int rounds = (int) log2(share_length); // Number of rounds for oblivious inequality
    // Step 1: Compute l_i = x_i XOR y_i XOR 1, 1 <= i < BSHARE_LENGTH
    BShare local_x_xor_y, remote_x_xor_y, lbs, rbs, shifted_lbs, shifted_rbs, l_diag, r_diag;
    for (int i=0; i<numElements; i++) {
        local_x_xor_y = c[i][leftcol] ^ c[i][rightcol];
        remote_x_xor_y = c[i][leftcol+1] ^ c[i][rightcol+1];
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
        l_diag |= (get_bit(c[i][rightcol], 0) ^ (BShare) 1);
        r_diag |= (get_bit(c[i][rightcol+1], 0) ^ (BShare) 1);
        // Local bit shares of the diagonal (store it in the second half of the array)
        local_bits[numElements+i] = and_b(l_diag, r_diag,
                                          c[i][leftcol], c[i][leftcol+1],
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
        result[i]=0;
        for (int j=0; j<share_length; j++)
            result[i] ^= get_bit(local_bits[i], j);
    }
    // Free memory
    free(local_bits); free(remote_bits);
}

// equi-join
void join_b(BShareTable input1, BShareTable input2, Predicate_B p, BShare result[]) {
  switch (p.operation) {
      case EQ:
          //printf("Equality predicate\n");
          join_eq_b(input1, input2, p.leftcol, p.rightcol, result);
          break;
      default:
          printf("Illegal operation. Only equality is supported for the moment.\n");
  }
}

// internal nested-loop equality join for boolean shares
PRIVATE void join_eq_b(BShareTable input1, BShareTable input2,
                        int leftcol, int rightcol, BShare result[]) {
  int i, j, k = 0;

  for (i = 0; i < input1.numRows; i++) {
    for (j = 0; j < input2.numRows; j++) {
      // generate rnums for the next equality
      result[k++] = eq_b(input1.content[i][leftcol], input1.content[i][leftcol + 1],
                         input2.content[j][rightcol], input2.content[j][rightcol + 1]);
    }
  }
}
 // batched join
void join_b_batch(BShareTable *input1, BShareTable *input2,
                    int start1, int end1, int start2, int end2,
                    Predicate_B p, BShare *remote, BShare *result) {

  switch (p.operation) {
      case EQ:
          //printf("Equality predicate\n");
          join_eq_b_batch(input1, input2,
                          start1, end1, start2, end2,
                          p.leftcol, p.rightcol, remote, result);
          break;
      case GEQ:
          //printf("Equality predicate\n");
          join_geq_b_batch(input1, input2,
                          start1, end1, start2, end2,
                          p.leftcol, p.rightcol, remote, result);
          break;
      default:
          printf("Illegal operation. Only equality is supported for the moment.\n");
  }

}

// batched version of internal nested-loop equality join for boolean shares
void join_eq_b_batch(BShareTable *input1, BShareTable *input2,
                        int start1, int end1, int start2, int end2,
                        int leftcol, int rightcol, BShare* remote,
                        BShare* result) {

  assert( (end1<=input1->numRows) && (end2<=input2->numRows) );

  int len1 = end1-start1;
  int len2 = end2-start2;
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);

  // compute bitwise x^y^1
  int k=0;
  for (int i=start1; i<end1; i++) {
    for (int j=start2; j<end2; j++) {
      result[k] = input1->content[i][leftcol] ^ input2->content[j][rightcol] ^ (~(BShare)0); // local share;
      remote[k] = input1->content[i][leftcol + 1] ^ input2->content[j][rightcol + 1] ^ (~(BShare)0); // remote share
      k++;
    }
  }

  // The result is stored in the (numbits/2) rightmost bits of result, res2 elements
  for (int l=0; l<numlevels; l++) {
    for (int i=0; i<len1*len2; i++) {
      result[i] = eq_b_level2(numbits >> l, result[i], remote[i]);
    }

    // Exchange results of logical and, except for the final round
    if (l != numlevels-1) {
      exchange_shares_array(result, remote, len1*len2);
    }
  }
}

// batched version of internal nested-loop equality join for boolean shares
void join_geq_b_batch(BShareTable *input1, BShareTable *input2,
                        int start1, int end1, int start2, int end2,
                        int leftcol, int rightcol, BShare* _remote,
                        BShare* result) {
    assert( (end1<=input1->numRows) && (end2<=input2->numRows) );
    BShare** c1 = input1->content;
    BShare** c2 = input2->content;

    int len1 = end1-start1;
    int len2 = end2-start2;

    int share_length=sizeof(BShare)*8; // The length of BShare in number of bits
    long numElements = len1 * len2;

    // Local and remote bits - We need two BShares per comparison
    unsigned long long *local_bits = (unsigned long long *) malloc(2 * numElements * sizeof(BShare));
    assert(local_bits!=NULL);
    unsigned long long *remote_bits = (unsigned long long *) malloc(2 * numElements * sizeof(BShare));
    assert(remote_bits!=NULL);
    int rounds = (int) log2(share_length); // Number of rounds for oblivious inequality
    // Step 1: Compute l_i = x_i XOR y_i XOR 1, 1 <= i < BSHARE_LENGTH
    BShare local_x_xor_y, remote_x_xor_y, lbs, rbs, shifted_lbs, shifted_rbs, l_diag, r_diag;
    int pos=0;
    for (int idx1=start1; idx1<end1; idx1++) {
        for (int idx2=start2; idx2<end2; idx2++) {
            local_x_xor_y = c2[idx2][rightcol] ^ c1[idx1][leftcol];
            remote_x_xor_y = c2[idx2][rightcol+1] ^ c1[idx1][leftcol+1];
            lbs = (local_x_xor_y ^ ~((BShare) 0));
            rbs = (remote_x_xor_y ^ ~((BShare) 0));
            // Compute first round
            shift_greater(lbs, rbs, 1, share_length, &shifted_lbs, &shifted_rbs);
            // Local bit shares of the 1st round
            local_bits[pos] = and_b(lbs, rbs,
                                    shifted_lbs, shifted_rbs,
                                    get_next_rb());
            // Set LSB to 'y_0 XOR 1' and compute the diagonal
            l_diag = unset_lsbs(local_x_xor_y, 1);
            r_diag = unset_lsbs(remote_x_xor_y, 1);
            l_diag |= (get_bit(c1[idx1][leftcol], 0) ^ (BShare) 1);
            r_diag |= (get_bit(c1[idx1][leftcol+1], 0) ^ (BShare) 1);
            // Local bit shares of the diagonal (store it in the second half of the array)
            local_bits[numElements+pos] = and_b(l_diag, r_diag,
                                                c2[idx2][rightcol], c2[idx2][rightcol+1],
                                                get_next_rb());
            pos++;
        }
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
        result[i]=0;
        for (int j=0; j<share_length; j++)
            result[i] ^= get_bit(local_bits[i], j);
        result[i] ^= (BitShare) 1;
    }
    // Free memory
    free(local_bits); free(remote_bits);
}

// Compares adjacent elements and populates 'distinct' array with secrets
// Expects an allocated 'distinct' array
void distinct(BShareTable* table, unsigned att_index, BitShare* distinct) {
  BShare** c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // First element of the sorted table is always in the set of distinct elements
  distinct[0] = eq_b(max, max, c[0][att_index], c[0][att_index+1])
                    ^ mask;
  // Compare elements in adjacent rows
  for (int i=0; i<table->numRows-1; i++) {
    distinct[i+1] = eq_b(c[i][att_index], c[i][att_index+1],
                         c[i+1][att_index], c[i+1][att_index+1])
                         ^ mask;
  }
}

// Same as distinct but works in batch mode
void distinct_batch(BShareTable* table, unsigned* key_indices, unsigned num_keys, 
                    BitShare* distinct, unsigned num_comparisons) {
  assert( (num_comparisons>0) && (num_comparisons<=table->numRows-1) );
  BShare **c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Allocate arrays for local and remote shares of equality results
  BShare *local = (BShare *) malloc(num_comparisons*num_keys*sizeof(BShare));
  assert(local!=NULL);
  BShare *remote = (BShare *) malloc(num_comparisons*num_keys*sizeof(BShare));
  assert(remote!=NULL);
  // First element of the sorted table is always in the set of distinct elements
  distinct[0] = get_rank() / 2;
  int num_bits = sizeof(BShare)*8;
  int num_levels = log2(num_bits);
  int total_comparisons = table->numRows - 1;
  // For each 'sliding' batch
  int next_start, batch_comp, offset, att_index;
  for (int i=0; i<total_comparisons; i+=num_comparisons) {
    // The start index of the next batch
    next_start = i + num_comparisons;
    // The number of comparisons in the current batch
    batch_comp = next_start <= total_comparisons ?
                               num_comparisons : total_comparisons-i;
    for (int j=0, k=i; j<batch_comp; j++, k++) {
      offset = j*num_keys;
      for (int idx=0; idx<num_keys; idx++) {
        att_index = key_indices[idx];
        // x_i ^ y_i ^ 1
        local[offset+idx] = c[k][att_index] ^ c[k+1][att_index] ^ max;
        remote[offset+idx] = c[k][att_index+1] ^ c[k+1][att_index+1] ^ max;
      }
    }
    // Apply equalities in bulk
    eq_bulk(num_levels, num_bits, local, remote, batch_comp*num_keys);
    exchange_shares_array(local, remote, batch_comp*num_keys);  // 1 round
    // AND all equality bits per pair of adjacent rows to compute the final bit bs
    and_b_all_group(local, remote, batch_comp, num_keys);   // log(num_keys) rounds
    // Set results
    for (int j=0, k=i+1; j<batch_comp; j++, k++) {
      distinct[k] = local[j] ^ mask;
    }
  }
  free(local);
  free(remote);
}

// Same as distinct_batch but returns BShare
PRIVATE void distinct_batch_incr(BShareTable* table, int start, int end,
                                    unsigned att_index, BShare* distinct) {
  int batch_size = end-start;
  int num_comparisons = start==0 ? batch_size-1 : batch_size;
  if (num_comparisons==0) { // We need one last comparison for the last element
    num_comparisons=1;
  }
  assert(num_comparisons>0);

  BShare **c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Allocate arrays for local and remote shares of equality results
  BShare *local = (BShare *) malloc(num_comparisons*sizeof(BShare));
  assert(local!=NULL);
  BShare *remote = (BShare *) malloc(num_comparisons*sizeof(BShare));
  assert(remote!=NULL);

  int num_bits = sizeof(BShare)*8;
  int num_levels = log2(num_bits);
  int pos;

  if (start==0) {
    // First element of the sorted table is always in the set of distinct
    distinct[0] = get_rank() / 2;
    pos = 1;
  }
  else {
    // We need one more comparison with the last element of the previous batch
    start -= 1;
    pos = 0;
  }

  // For each 'sliding' batch
  for (int i=start; i<end-1; i+=num_comparisons) {
    for (int j=0, k=i; j<num_comparisons; j++, k++) {
      // x_i ^ y_i ^ 1
      local[j] = c[k][att_index] ^ c[k+1][att_index] ^ max;
      remote[j] = c[k][att_index+1] ^ c[k+1][att_index+1] ^ max;
    }
    // Apply equalities in bulk
    eq_bulk(num_levels, num_bits, local, remote, num_comparisons);
    // Set results
    for (int j=pos, k=0; j<batch_size; j++, k++) {
      distinct[j] = local[k] ^ mask;
    }
  }
  // Free memory
  free(local); free(remote);
}

// Computes distinct with a linear scan taking a previously computed
// selected bit into account
void distinct_linear(BShareTable* table, unsigned* key_indices, int num_keys,
                    BShare* selected_b) {
  BShare** c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Fetch the second boolean share of each 'selected' bit -- 1 round
  BShare *remote_selected_b = (BShare *) malloc((table->numRows)*sizeof(BShare));
  assert(remote_selected_b!=NULL);
  exchange_shares_array(selected_b, remote_selected_b, table->numRows);

  int len = table->numCols/2 + 1;   // Number of attributes plus propagated 
                                    // 'selected' bit
  int rank = get_rank();
  int succ_rank = get_succ();
  BShare bs1, bs2;
  BShare local[len], remote[len];
  BShare local_bits[num_keys], remote_bits[num_keys];

  for (int i=0, k=0; i<table->numRows-1; i++, k+=2) {
    // Compute a single bit bs that denotes whether the adjacent rows c[i],
    // c[i+1] have the same id
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
    // bs denotes whether the adjacent rows c[i], c[i+1] have the same id
    bs1 = local_bits[0];
    bs2 = remote_bits[0];
    // Compute bn = (bs OR NOT selected_b[i]) = NOT(NOT bs AND selected_b[i])
    BShare bn1 = bs1 ^ mask,
           bn2 = bs2 ^ mask;
    bn1 = and_b(bn1, bn2, selected_b[i] & mask, remote_selected_b[i] & mask,
                get_next_rb()) ^ mask;
    bn2 = exchange_shares(bn1);         // 1 round
    bn1 &= mask;                        // Keep LSB only
    bn2 &= mask;                        // Keep LSB only
    BShare b1 = -bn1;                   // Set all bits equal to LSB of bn1
    BShare b2 = -bn2;                   // Set all bits equal to LSB of bn2
    // Compute new_c[i] = b * dummy_row + (1-b) * row1
    for (int j=0; j<table->numCols; j+=2) {
      local[j/2] = and_b(b1, b2, max, max, get_next_rb());
      local[j/2] ^= and_b(~b1, ~b2, c[i][j], c[i][j+1], get_next_rb());
    }
    // Compute cond = bs AND selected_b[i] to propagate 'selected' bit
    BShare cond1 = and_b(bs1, bs2, selected_b[i], remote_selected_b[i],
                         get_next_rb());
    BShare cond2 = exchange_shares(cond1);  // 1 round
    cond1 = -cond1;                         // Set all bits equal to LSB
    cond2 = -cond2;                         // Set all bits equal to LSB
    // Compute cond * selected_b[i] + (1-cond)*selected_b[i+1]
    local[len-1] = and_b(cond1, cond2, selected_b[i], remote_selected_b[i],
                         get_next_rb());
    local[len-1] ^= and_b(~cond1, ~cond2,
                          selected_b[i+1], remote_selected_b[i+1],
                          get_next_rb());
    // Fetch remote boolean shares
    exchange_shares_array(local, remote, len);    // 1 round
    // Set c[i] = new_c[i]
    for (int j=0; j<table->numCols; j+=2) {
      c[i][j] = local[j/2];
      c[i][j+1] = remote[j/2];
    }
    // Propagate 'selected' bit
    selected_b[i+1] = local[len-1];
    remote_selected_b[i+1] = remote[len-1];

    // set selected_bit[i] = 0 if propagated
    // Compute ~cond * selected_b[i]
    selected_b[i] = and_b(~cond1, ~cond2, selected_b[i], remote_selected_b[i],
                         get_next_rb());
    exchange_shares_array(selected_b, remote_selected_b, table->numRows);
  }

  // Make sure we mask the last row if it's not selected
  // 1. Compute composite bit bl = NOT(bs) AND NOT(selected_b)
  BShare bl1 = and_b(bs1 ^ mask, bs2 ^ mask,
                     selected_b[table->numRows-1] ^ mask,
                     remote_selected_b[table->numRows-1] ^ mask,
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
  // Final shuffle iff not followed by ORDER_BY
  free(remote_selected_b);
}

// Applies geq() to adjacent rows in batch
void adjacent_geq(BShareTable* table, unsigned att_index1, unsigned att_index2,
                  BitShare* result, unsigned num_comparisons, int swap) {
    assert( (num_comparisons>0) && (num_comparisons<=table->numRows-1) );
    BShare **c = table->content;
    // Local and remote bits - We need two BShares per comparison
    unsigned long long *local_bits = (unsigned long long *) malloc(2 * num_comparisons * sizeof(BShare));
    assert(local_bits!=NULL);
    unsigned long long *remote_bits = (unsigned long long *) malloc(2 * num_comparisons * sizeof(BShare));
    assert(remote_bits!=NULL);
    int share_length = sizeof(BShare)*8; // The length of BShare in number of bits
    int rounds = (int) log2(share_length); // Number of rounds for oblivious inequality
    int total_comparisons = table->numRows - 1, batch_comp, next_start;
    for (int p=0; p<total_comparisons; p+=num_comparisons) {
        // The start index of the next batch
        next_start = p + num_comparisons;
        // The number of comparisons in the current batch
        batch_comp = next_start <= total_comparisons ?
                     num_comparisons : total_comparisons-p;
        // Step 1: Compute l_i = x_i XOR y_i XOR 1, 1 <= i < BSHARE_LENGTH
        BShare local_x_xor_y, remote_x_xor_y, lbs, rbs, shifted_lbs, shifted_rbs, l_diag, r_diag;
        for (int i=0, k=p; i<batch_comp; i++, k++) {
            local_x_xor_y = c[k+1-swap][att_index2] ^ c[k+swap][att_index1];
            remote_x_xor_y = c[k+1-swap][att_index2+1] ^ c[k+swap][att_index1+1];
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
            l_diag |= (get_bit(c[k+swap][att_index1], 0) ^ (BShare) 1);
            r_diag |= (get_bit(c[k+swap][att_index1+1], 0) ^ (BShare) 1);
            // Local bit shares of the diagonal (store it in the second half of the array)
            local_bits[batch_comp+i] = and_b(l_diag, r_diag,
                                             c[k+1-swap][att_index2], c[k+1-swap][att_index2+1],
                                             get_next_rb());
        }
        // Fetch remote bit shares
        exchange_shares_array_u(local_bits, remote_bits, batch_comp*2);
        // Compute the remaining rounds
        for (int r=2; r<=rounds; r++) {
            for (int i=0; i<batch_comp; i++) {
                lbs = local_bits[i];
                rbs = remote_bits[i];
                shift_greater(lbs, rbs, r, share_length, &shifted_lbs, &shifted_rbs);
                local_bits[i] = and_b(lbs, rbs, shifted_lbs, shifted_rbs, get_next_rb());
            }
            // Fetch remote shares
            exchange_shares_array_u(local_bits, remote_bits, batch_comp);
        }
        // Do a final AND with the diagonal
        for (int i=0; i<batch_comp; i++) {
            lbs = local_bits[i];
            rbs = remote_bits[i];
            l_diag = local_bits[batch_comp+i];
            r_diag = remote_bits[batch_comp+i];
            local_bits[i] = and_b(set_bit(lbs>>1, share_length-1), set_bit(rbs>>1, share_length-1),
                                  l_diag, r_diag,
                                  get_next_rb());
        }
        // Fetch remote shares
        exchange_shares_array_u(local_bits, remote_bits, batch_comp);
        // XOR all bits of each result and negate it
        for (int i=0, k=p; i<batch_comp; i++, k++) {
            result[k]=0;
            for (int j=0; j<share_length; j++)
                result[k] ^= get_bit(local_bits[i], j);
            result[k] ^= (BitShare) 1;
        }
    }
    // Free memory
    free(local_bits); free(remote_bits);
}

// Oblivious semi-join -- relies on join_eq_b_batch()
void in(BShareTable* left, BShareTable* right, unsigned left_index,
        unsigned right_index, BShare* in, unsigned num_rows_left) {
  // Right input must be a power of two
  assert( ceil(log2(right->numRows)) == floor(log2(right->numRows)));
  BShare mask=1;
  int num_levels = log2(right->numRows);
  // Total number of equalities per batch
  long num_comparisons = num_rows_left * right->numRows;
  // Allocate arrays for local and remote shares of equality results
  BShare *local = (BShare *) malloc(num_comparisons*sizeof(BShare));
  assert(local!=NULL);
  BShare *remote = (BShare *) malloc(num_comparisons*sizeof(BShare));
  assert(remote!=NULL);
  int next, left_end, num_comp_batch;
  // For all batches
  for (int i=0; i<left->numRows; i+=num_rows_left) {
    next = i + num_rows_left;
    left_end = next <= left->numRows ? next : left->numRows;
    num_comp_batch = (left_end - i) * right->numRows; // Comparisons in batch
    // 1. Do the first round of the semi-join for all 'num_rows_left' in bulk
    join_eq_b_batch(left, right, i, left_end, 0, right->numRows,
                    left_index, right_index, remote, local);
    exchange_shares_array(local, remote, num_comp_batch);
    // 2. Do the second round of the semi-join for all 'num_rows_left' in bulk
    for (int j=0; j<num_comp_batch; j++){
      // Compute eq_result ^ 1
      local[j] ^= mask;
      local[j] &= mask;
      remote[j] ^= mask;
      remote[j] &= mask;
    }
    // Compute equalities in bulk
    for (int l=0; l<num_levels; l++) {
      // For each row on the left
      for (int j=0; j<num_comp_batch; j+=right->numRows) {
        in_level(right->numRows >> l, &local[j], &remote[j]);
      }
      // Exchange results of logical AND, except for the final round
      if (l != num_levels-1) {
        exchange_shares_array(local, remote, num_comp_batch);
      }
    }
    // Set IN results
    for (int k=i, j=0; k<left_end; k++, j+=right->numRows) {
      in[k] = (local[j] ^ mask);
    }
  }
  free(local);
  free(remote);
}

// operates like in() but also takes into account a selected bit
// (sel_index) in the right table
void in_sel_right(BShareTable* left, BShareTable* right,
        unsigned left_index, unsigned right_index, unsigned sel_index,
        BShare* res_in, unsigned num_rows_left) {

  // Right input must be a power of two
  assert( ceil(log2(right->numRows)) == floor(log2(right->numRows)));
  BShare mask=1;
  int num_levels = log2(right->numRows);
  // Total number of equalities per batch
  long num_comparisons = num_rows_left * right->numRows;
  // Allocate arrays for local and remote shares of equality results
  BShare *local = (BShare *) malloc(num_comparisons*sizeof(BShare));
  assert(local!=NULL);
  BShare *remote = (BShare *) malloc(num_comparisons*sizeof(BShare));
  assert(remote!=NULL);
  int next, left_end, num_comp_batch;

  // For all batches
  for (int i=0; i<left->numRows; i+=num_rows_left) {
    next = i + num_rows_left;
    left_end = next <= left->numRows ? next : left->numRows;
    num_comp_batch = (left_end - i) * right->numRows; // Comparisons in batch

    // 1. Do the first round of the semi-join for all 'num_rows_left' in bulk
    join_eq_b_batch(left, right, i, left_end, 0, right->numRows,
                    left_index, right_index, remote, local);
    exchange_shares_array(local, remote, num_comp_batch);

    // Take into account the selected bit
    // We need to compute a logical and of local with selected
    for (long i=0; i<num_comp_batch; i++) {
      local[i] = and_b(local[i], remote[i],
                      right->content[i % right->numRows][sel_index],
                      right->content[i % right->numRows][sel_index + 1],
                      get_next_rb());
    }
    exchange_shares_array(local, remote, num_comp_batch);

    // 2. Do the second round of the semi-join for all 'num_rows_left' in bulk
    for (int j=0; j<num_comp_batch; j++){
      // Compute eq_result ^ 1
      local[j] ^= mask;
      local[j] &= mask;
      remote[j] ^= mask;
      remote[j] &= mask;
    }
    // Compute equalities in bulk
    for (int l=0; l<num_levels; l++) {
      // For each row on the left
      for (int j=0; j<num_comp_batch; j+=right->numRows) {
        in_level(right->numRows >> l, &local[j], &remote[j]);
      }
      // Exchange results of logical AND, except for the final round
      if (l != num_levels-1) {
        exchange_shares_array(local, remote, num_comp_batch);
      }
    }
    // Set IN results
    for (int k=i, j=0; k<left_end; k++, j+=right->numRows) {
      res_in[k] = (local[j] ^ mask);
    }
  }
  free(local);
  free(remote);


}

// Used by in()
PRIVATE void in_level(int num_elements, BShare* z1, BShare* z2) {
  const BShare mask = 1;
  for (int i=0, j=0; i<num_elements; i+=2, j++) {
    BShare bx1 = z1[i];
    BShare bx2 = z2[i];
    BShare by1 = z1[i+1];
    BShare by2 = z2[i+1];
    // store the result (and's LSB) in the result's jth bit
    BShare out = and_b(bx1, bx2, by1, by2, get_next_rb());
    z1[j] = (out & mask);
  }
}

// Used by in() and distinct_batch()
void eq_bulk(int num_levels, int num_bits, BShare* local,
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

// Groups rows and counts the number of rows per group using an odd-even merge
// Takes into account the selected rows from previous operations
// The given BShareTable must be sorted on the group-by key(s)
void group_by_count_sel_odd_even(BShareTable* table, unsigned* key_indices, int num_keys,
                        unsigned batch_size, BShare* selected_b, AShare* selected, 
                        BShare* rb, AShare* ra) {
  // Make sure the number of rows is a power of two
  assert(ceil(log2(table->numRows)) == floor(log2(table->numRows)));
  // Make sure batch_size is at most equal to the input size
  batch_size = batch_size > table->numRows ? table->numRows : batch_size;
  BShare** c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Fetch the second arithmetic share of each 'selected' bit (1 round)
  AShare *remote_selected = (AShare *) malloc((table->numRows)*sizeof(AShare));
  assert(remote_selected!=NULL);
  exchange_a_shares_array(selected, remote_selected, table->numRows);
  // Fetch the second boolean share of each 'selected' bit (1 round)
  BShare *remote_selected_b = (BShare *) malloc((table->numRows)*sizeof(BShare));
  assert(remote_selected_b!=NULL);
  exchange_shares_array(selected_b, remote_selected_b, table->numRows);
  // Allocate memory for local and remote shares
  BShare  *local_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)),  // local shares of equality bits
          *remote_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)); // remote shares of equality bits
  assert(local_bits!=NULL); assert(remote_bits!=NULL);
  AShare  *ab1= (AShare* ) malloc(batch_size*sizeof(AShare)),  // local shares of arithmetic equality bits
          *ab2= (AShare* ) malloc(batch_size*sizeof(AShare));  // remote shares of arithmetic equality bits
  assert(ab1!=NULL); assert(ab2!=NULL);
  int len = table->numCols/2 + 1;  // Number of attributes plus new count
  BShare  *local= (BShare  *) malloc(batch_size*len*sizeof(BShare)),  // local shares of multiplexed rows
          *remote= (BShare  *) malloc(batch_size*len*sizeof(BShare)); // remote shares of multiplexed rows
  assert(local!=NULL); assert(remote!=NULL);
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int rank = get_rank();
  int succ_rank = get_succ();
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
      // Compute new_c[i+dist] = bn * dummy_row + (1-bn) * c[i+dist] to multiplex rows
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
      // Compute arithmetic shares of the equality bits using the boolean shares
      convert_single_bit_array(local_bits, &ra[num_comparisons], &rb[num_comparisons],
                               num_pairs, ab1);   // 2 rounds
      exchange_shares_array(ab1, ab2, num_pairs); // 1 round
      for (int i=start, k=0; i<end; i++, k++) {
        AShare a_bs1 = ab1[k], a_bs2 = ab2[k];
        // Compute new_cnt = bn*(selected[i] + selected[i+dist]) + (1-bn)*selected[i]
        AShare local_count = selected[i] + selected[i+dist],
               remote_count = remote_selected[i] + remote_selected[i+dist];
        offset = (k+1) * len;
        local[offset-1] = mul(a_bs1, a_bs2, local_count, remote_count,
                              get_next_r());
        local[offset-1] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, selected[i],
                               remote_selected[i], get_next_r());
      }
      // Fetch remote boolean and arithmetic shares in a single round
      // NOTE: This works because BShare and AShare are both of the same type
      exchange_shares_array(local, remote, num_pairs*len);    // 1 round
      // Update rows in place
      for (int i=start, k=0; i<end; i++, k++) {
        offset = k*len;
        // Set c[i+dist] = new_c[i+dist]
        for (int j=0; j<table->numCols; j+=2) {
          c[i+dist][j] = local[offset+j/2];
          c[i+dist][j+1] = remote[offset+j/2];
        }

        // Set selected[i] = new_cnt
        selected[i] = local[offset+len-1];
        remote_selected[i] = remote[offset+len-1];

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
    int rindex = num_comparisons;  // The index of the first available random number
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
        convert_single_bit_array(local, &ra[rindex], &rb[rindex],
                                 comp, ab1);   // 2 rounds
        exchange_shares_array(ab1, ab2, comp); // 1 round
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
            // Compute selected[i] = b*max + (1-b)*selected[i]
            AShare a_b1 = ab1[k], a_b2 = ab2[k];
            local[k+offset] = mul(a_b1, a_b2, max, max, get_next_r());
            local[k+offset] += mul(rank%2 - a_b1, succ_rank%2 - a_b2,
                                   selected[i], remote_selected[i], get_next_r());
        }
        exchange_shares_array(local, remote, comp+offset); // 1 round
        rindex += comp;
        if (start==0) {  // Update first row
            for (int j=0; j<table->numCols; j+=2) {
                c[0][j] = local[j/2];
                c[0][j+1] = remote[j/2];
            }
        }
        // Update counts
        for (int i=start, k=0; i<end; i++, k++) {
            selected[i] = local[k+offset];
            remote_selected[i] = remote[k+offset];
        }
    }
    // NOTE: Final shuffle iff not followed by ORDER_BY
    free(remote_selected); free(remote_selected_b);
    free(ab1); free(ab2);
    free(local_bits); free(remote_bits);
    free(local); free(remote);
}

// Groups rows and counts the number of rows per group using an odd-even merge
// The given BShareTable must be sorted on the group-by key(s)
void group_by_count_odd_even(BShareTable* table, unsigned* key_indices, int num_keys,
                             unsigned batch_size, AShare* counters,
                             AShare* remote_counters, BShare* rb, AShare* ra) {
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
  AShare  *ab1= (AShare  *) malloc(batch_size*sizeof(AShare)),  // local shares of arithmetic equality bits
          *ab2= (AShare  *) malloc(batch_size*sizeof(AShare));  // remote shares of arithmetic equality bits
  assert(ab1!=NULL); assert(ab2!=NULL);
  int len = table->numCols/2 + 1;  // Number of attributes plus new count
  BShare  *local= (BShare  *) malloc(batch_size*len*sizeof(BShare)),  // local shares of multiplexed rows
          *remote= (BShare  *) malloc(batch_size*len*sizeof(BShare)); // remote shares of multiplexed rows
  assert(local!=NULL); assert(remote!=NULL);
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int rank = get_rank();
  int succ_rank = get_succ();
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
      // Compute arithmetic shares of the equality bits using the boolean shares
      convert_single_bit_array(local_bits, &ra[num_comparisons], &rb[num_comparisons], 
                               num_pairs, ab1);   // 2 rounds
      exchange_shares_array(ab1, ab2, num_pairs); // 1 round
      for (int i=start, k=0; i<end; i++, k++) {
        AShare a_bs1 = ab1[k], a_bs2 = ab2[k];
        // Compute new_cnt = bs*(counters[i] + counters[i+dist]) + (1-bs)*counters[i]
        AShare local_count = counters[i] + counters[i+dist],
               remote_count = remote_counters[i] + remote_counters[i+dist];
        offset = (k+1) * len;
        local[offset-1] = mul(a_bs1, a_bs2, local_count, remote_count,
                              get_next_r());
        local[offset-1] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, counters[i],
                               remote_counters[i], get_next_r());
      }
      // Fetch remote boolean and arithmetic shares in a single round
      // NOTE: This works because BShare and AShare are both of the same type
      exchange_shares_array(local, remote, num_pairs*len);    // 1 round
      // Update rows in place
      for (int i=start, k=0; i<end; i++, k++) {
        offset = k*len;
        // Set c[i+dist] = new_c[i+dist]
        for (int j=0; j<table->numCols; j+=2) {
          c[i+dist][j] = local[offset+j/2];
          c[i+dist][j+1] = remote[offset+j/2];
        }
        // Set selected[i] = new_cnt
        counters[i] = local[offset+len-1];
        remote_counters[i] = remote[offset+len-1];
      }
      num_comparisons += num_pairs;
      // Update batch boundaries
      start = end;
      step = end + batch_size;
      end = (step <= total_comparisons ? step : total_comparisons);
      num_pairs = end - start;
    } while (start<end);
  }
  // Do a final pass to mask counts of masked rows
  int rindex = num_comparisons;  // The index of the first available random number
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
    convert_single_bit_array(local_bits, &ra[rindex], &rb[rindex],
                             comp, ab1);   // 2 rounds
    exchange_shares_array(ab1, ab2, comp); // 1 round
    // Compute new_cnt[start] = b*max + (1-b)*selected[start]
    for (int i=start, k=0; i<end; i++, k++) {
      AShare a_b1 = ab1[k], a_b2 = ab2[k];
      local_bits[k] = mul(a_b1, a_b2, max, max, get_next_r());
      local_bits[k] += mul(rank%2 - a_b1, succ_rank%2 - a_b2,
                           counters[i], remote_counters[i], get_next_r());
    }
    exchange_shares_array(local_bits, remote_bits, comp); // 1 round
    rindex += comp;
    // Update counts
    for (int i=start, k=0; i<end; i++, k++) {
      counters[i] = local_bits[k];
      remote_counters[i] = remote_bits[k];
    }
  }
  // NOTE: Final shuffle iff not followed by ORDER_BY
  free(ab1); free(ab2);
  free(local_bits); free(remote_bits);
  free(local); free(remote);
}

// Groups rows and counts the number of rows per group
// Takes into account the selected rows from previous operations
// The given BShareTable must be sorted on the group-by key(s)
void group_by_count(BShareTable* table, unsigned* key_indices, int num_keys,
                    BShare* selected_b, AShare* selected, BShare* rb,
                    AShare* ra) {
  BShare** c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Fetch the second arithmetic share of each 'selected' bit -- 1 round
  AShare *remote_selected = (AShare *) malloc((table->numRows)*sizeof(AShare));
  assert(remote_selected!=NULL);
  exchange_a_shares_array(selected, remote_selected, table->numRows);
  // Fetch the second boolean share of each 'selected' bit -- 1 round
  BShare *remote_selected_b = (BShare *) malloc((table->numRows)*sizeof(BShare));
  assert(remote_selected_b!=NULL);
  exchange_shares_array(selected_b, remote_selected_b, table->numRows);
  // Scan table and update counts by adding 'selected' bits
  AShare ab1[2], ab2[2];
  BShare bb[2];
  int len = table->numCols/2 + 3;   // Number of attributes plus old and new
                                    // counts plus propagated 'selected' bit
  int rank = get_rank();
  int succ_rank = get_succ();
  BShare bs1, bs2;
  BShare local[len], remote[len];
  BShare local_bits[num_keys], remote_bits[num_keys];
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
    BShare bn1 = bs1 ^ mask,
           bn2 = bs2 ^ mask;
    bn1 = and_b(bn1, bn2, selected_b[i] & mask, remote_selected_b[i] & mask,
                get_next_rb()) ^ mask;
    bn2 = exchange_shares(bn1);         // 1 round
    bn1 &= mask;                        // Keep LSB only
    bn2 &= mask;                        // Keep LSB only
    BShare b1 = -bn1;                   // Set all bits equal to LSB of bn1
    BShare b2 = -bn2;                   // Set all bits equal to LSB of bn2
    // Compute new_c[i] = b * dummy_row + (1-b) * row1
    for (int j=0; j<table->numCols; j+=2) {
      local[j/2] = and_b(b1, b2, max, max, get_next_rb());
      local[j/2] ^= and_b(~b1, ~b2, c[i][j], c[i][j+1], get_next_rb());
    }
    // Compute arithmetic shares from boolean shares
    bb[0] = bs1;
    bb[1] = bn1;
    convert_single_bit_array(bb, &ra[k], &rb[k], 2, ab1);   // 1 round
    exchange_shares_array(ab1, ab2, 2);                     // 1 round
    AShare a_bs1 = ab1[0], a_bs2 = ab2[0];
    AShare a_bn1 = ab1[1], a_bn2 = ab2[1];
    // Compute new_cnt = bs*(selected[i] + selected[i+1]) + (1-bs)*selected[i+1]
    AShare local_count = selected[i] + selected[i+1],
           remote_count = remote_selected[i] + remote_selected[i+1];
    local[len-3] = mul(a_bs1, a_bs2, local_count, remote_count,
                       get_next_r());
    local[len-3] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, selected[i+1],
                        remote_selected[i+1], get_next_r());
    // Compute bn * max + (1-bn) * selected[i] to allow masking 'selected' bit
    local[len-2] = mul(a_bn1, a_bn2, max, max, get_next_r());
    local[len-2] += mul(rank%2 - a_bn1, succ_rank%2 - a_bn2,
                        selected[i], remote_selected[i], get_next_r());
    // Compute cond = bs AND selected_b[i] to propagate 'selected' bit
    BShare cond1 = and_b(bs1, bs2, selected_b[i], remote_selected_b[i],
                         get_next_rb());
    BShare cond2 = exchange_shares(cond1);  // 1 round
    cond1 = -cond1;                         // Set all bits equal to LSB
    cond2 = -cond2;                         // Set all bits equal to LSB
    // Compute cond * selected_b[i] + (1-cond)*selected_b[i+1]
    local[len-1] = and_b(cond1, cond2, selected_b[i], remote_selected_b[i],
                         get_next_rb());
    local[len-1] ^= and_b(~cond1, ~cond2,
                          selected_b[i+1], remote_selected_b[i+1],
                          get_next_rb());
    // Fetch remote boolean and arithmetic shares
    // NOTE: This works because BShare and AShare are both long long
    exchange_shares_array(local, remote, len);    // 1 round
    // Set c[i] = new_c[i]
    for (int j=0; j<table->numCols; j+=2) {
      c[i][j] = local[j/2];
      c[i][j+1] = remote[j/2];
    }
    // Propagate 'selected' bit
    selected_b[i+1] = local[len-1];
    remote_selected_b[i+1] = remote[len-1];
    // Update 'selected' bit
    selected[i] = local[len-2];
    remote_selected[i] = remote[len-2];
    // Set selected[i+1] = new_cnt
    selected[i+1] = local[len-3];
    remote_selected[i+1] = remote[len-3];
  }
  // Make sure we mask the last row if it's not selected
  // 1. Compute composite bit bl = NOT(bs) AND NOT(selected_b)
  BShare bl1 = and_b(bs1 ^ mask, bs2 ^ mask,
                     selected_b[table->numRows-1] ^ mask,
                     remote_selected_b[table->numRows-1] ^ mask,
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
  // Final shuffle iff not followed by ORDER_BY
  free(remote_selected); free(remote_selected_b);
}

// Groups rows and computes min and max of the specified attributes per group
// Takes into account the selected rows from previous operations
// The given BShareTable must be sorted on the group-by key(s)
void group_by_min_max_sel(BShareTable* table, BShare* selected,
                          unsigned min_att, unsigned max_att,
                          unsigned* key_indices, int num_keys) {
  BShare** c = table->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Fetch the second arithmetic share of each 'selected' bit -- 1 round
  AShare *remote_selected = (AShare *) malloc((table->numRows)*sizeof(AShare));
  assert(remote_selected!=NULL);
  exchange_a_shares_array(selected, remote_selected, table->numRows);
  int len = table->numCols/2 + 4;   // Number of attributes plus new min and max
                                    // plus masked and propagated selected bits
  BShare bs1, bs2, bb1[2], bb2[2];  // Bit shares
  BShare local[len], remote[len];
  BShare local_bits[num_keys], remote_bits[num_keys];
  // Scan table and update min/max
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
    // Compute bn = (bs OR NOT selected[i]) = NOT(NOT bs AND selected[i])
    // bn is used to mask the i-th row when "it is not in the same group with
    // the next row or it is not selected"
    BShare bn1 = bs1 ^ mask,
           bn2 = bs2 ^ mask;
    bn1 = and_b(bn1, bn2, selected[i] & mask, remote_selected[i] & mask,
                get_next_rb()) ^ mask;
    bn2 = exchange_shares(bn1);         // 1 round
    bn1 &= mask;                        // Keep LSB only
    bn2 &= mask;                        // Keep LSB only
    BShare b1 = -bn1;                   // Set all bits equal to LSB of bn1
    BShare b2 = -bn2;                   // Set all bits equal to LSB of bn2
    // Compute new_c[i] = b * dummy_row + (1-b) * row1
    for (int j=0; j<table->numCols; j+=2) {
      local[j/2] = and_b(b1, b2, max, max, get_next_rb());
      local[j/2] ^= and_b(~b1, ~b2, c[i][j], c[i][j+1], get_next_rb());
    }
    // Compute bmin = min_att[i] ?< min_att[i+1]
    bb1[0] = greater(c[i+1][min_att], c[i+1][min_att+1], c[i][min_att],
                     c[i][min_att+1]);
    // Compute bmax = max_att[i] ?> max_att[i+1]
    bb1[1] = greater(c[i][max_att], c[i][max_att+1], c[i+1][max_att],
                     c[i+1][max_att+1]);
    exchange_shares_array(bb1, bb2, 2);     // 1 round
    BShare bmin1 = bb1[0];
    BShare bmin2 = bb2[0];
    BShare bmax1 = bb1[1];
    BShare bmax2 = bb2[1];
    // Compute bl = selected[i] AND selected[i+1] ("both selected")
    bb1[0] = and_b(selected[i], remote_selected[i], selected[i+1],
                   remote_selected[i+1], get_next_rb())
                 & mask;
    // Compute bf = selected[i] AND NOT selected[i+1] ("first selected")
    bb1[1] = and_b(selected[i], remote_selected[i], selected[i+1] ^ mask,
                   remote_selected[i+1] ^ mask, get_next_rb())
                 & mask;
    exchange_shares_array(bb1, bb2, 2);         // 1 round
    BShare bl1 = bb1[0], bl2 = bb2[0], bf1 = bb1[1], bf2 = bb2[1];
    // Compute bc = bs AND bl ("same group and both selected")
    bb1[0] = and_b(bs1, bs2, bl1, bl2, get_next_rb()) & mask;
    // Compute br = bs AND bf ("same group and first selected")
    bb1[1] = and_b(bs1, bs2, bf1, bf2, get_next_rb()) & mask;
    exchange_shares_array(bb1, bb2, 2);       // 1 round
    BShare bc1 = bb1[0], bc2 = bb2[0], br1 = bb1[1], br2 = bb2[1];;
    // Compute composite bits:
    //  - bcx = bc AND bmax
    //  - bcn = bc AND bmin
    bb1[0] = and_b(bc1, bc2, bmax1, bmax2, get_next_rb()) & mask;
    bb1[1] = and_b(bc1, bc2, bmin1, bmin2, get_next_rb()) & mask;
    exchange_shares_array(bb1, bb2, 2); // 1 round
    BShare bcx1 = bb1[0], bcx2 = bb2[0], bcn1 = bb1[1], bcn2 = bb2[1];
    // Compute composite bits:
    //  - bu = bcx OR br = NOT (NOT bcx AND NOT br)
    //  - bp = bcn OR br = NOT (NOT bcn AND NOT br)
    bb1[0] = (and_b(bcx1 ^ mask, bcx2 ^ mask, br1 ^ mask, br2 ^ mask,
                    get_next_rb()) ^ mask)
                 & mask;
    bb1[1] = (and_b(bcn1 ^ mask, bcn2 ^ mask, br1 ^ mask, br2 ^ mask,
                    get_next_rb()) ^ mask)
                 & mask;
    exchange_shares_array(bb1, bb2, 2); // 1 round
    BShare bu1 = bb1[0], bu2 = bb2[0], bp1 = bb1[1], bp2 = bb2[1];
    // Compute new_max = bu*c[i] + (1-bu)*c[i+1]
    local[len-4] = and_b(-bu1, -bu2, c[i][max_att], c[i][max_att+1],
                         get_next_rb());
    local[len-4] ^= and_b(~(-bu1), ~(-bu2), c[i+1][max_att], c[i+1][max_att+1],
                         get_next_rb());
    // Compute new_min = bp*c[i] + (1-bp)*c[i+1]
    local[len-3] = and_b(-bp1, -bp2, c[i][min_att], c[i][min_att+1],
                         get_next_rb());
    local[len-3] ^= and_b(~(-bp1), ~(-bp2), c[i+1][min_att], c[i+1][min_att+1],
                         get_next_rb());
    // Compute bn * max + (1-bn) * selected[i] to mask the selected bit when
    // the "i-th row is in the same group with the next or it is not selected"
    local[len-2] = and_b(b1, b2, max, max, get_next_rb());
    local[len-2] ^= and_b(~b1, ~b2, selected[i], remote_selected[i],
                          get_next_rb());
    // Compute cond = bs AND selected[i] to propagate 'selected' bit when the
    // "i-th row is selected and belongs to the same group with the next row"
    BShare cond1 = and_b(bs1, bs2, selected[i], remote_selected[i],
                         get_next_rb());
    BShare cond2 = exchange_shares(cond1);  // 1 round
    cond1 = -cond1;                         // Set all bits equal to LSB
    cond2 = -cond2;                         // Set all bits equal to LSB
    // Compute cond * selected[i] + (1-cond)*selected[i+1]
    local[len-1] = and_b(cond1, cond2, selected[i], remote_selected[i],
                         get_next_rb());
    local[len-1] ^= and_b(~cond1, ~cond2, selected[i+1], remote_selected[i+1],
                          get_next_rb())
                        & mask;
    // Fetch remote shares
    exchange_shares_array(local, remote, len);    // 1 round
    // Set c[i] = new_c[i]
    for (int j=0; j<table->numCols; j+=2) {
      c[i][j] = local[j/2];
      c[i][j+1] = remote[j/2];
    }
    // Propagate 'selected' bit
    selected[i+1] = local[len-1];
    remote_selected[i+1] = remote[len-1];
    // Update previous 'selected' bit
    selected[i] = local[len-2];
    remote_selected[i] = remote[len-2];
    // Set new_min
    c[i+1][min_att] = local[len-3];
    c[i+1][min_att+1] = remote[len-3];
    // Set new_max
    c[i+1][max_att] = local[len-4];
    c[i+1][max_att+1] = remote[len-4];
  }
  // Make sure we mask the last row if it's not selected
  // 1. Compute composite bit bk = NOT(bs) AND NOT(selected)
  BShare bk1 = and_b(bs1 ^ mask, bs2 ^ mask,
                     selected[table->numRows-1] ^ mask,
                     remote_selected[table->numRows-1] ^ mask,
                     get_next_rb())
                   & mask;
  BShare bk2 = exchange_shares(bk1);            // 1 round
  BShare b1 = -bk1;
  BShare b2 = -bk2;
  // 2. Compute new_c[i] = bk * dummy_row + (1-bk) * last_row
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

// Same as group_by_min_max_sel() but uses odd-even aggregation
void group_by_min_max_sel_odd_even(BShareTable* table, unsigned batch_size, 
                                   BShare* selected_b,
                                   unsigned min_att, unsigned max_att,
                                   unsigned* key_indices, int num_keys) {
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
  BShare  *local_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)),  // local shares of equality bits
          *remote_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)); // remote shares of equality bits
  assert(local_bits!=NULL); assert(remote_bits!=NULL);
  BShare  *bits= (BShare  *) malloc(batch_size*sizeof(BShare)),  // local shares of intermediate bits
          *rbits= (BShare  *) malloc(batch_size*sizeof(BShare)); // remote shares of intermediate bits
  assert(bits!=NULL); assert(rbits!=NULL);
  BShare  *bcns= (BShare *) malloc(batch_size*sizeof(BShare)),  // local shares of intermediate bits
          *rbcns= (BShare *)malloc(batch_size*sizeof(BShare)); // remote shares of intermediate bits
  assert(bcns!=NULL); assert(rbcns!=NULL);
  BitShare  *bmins= (BitShare  *) malloc(batch_size*sizeof(BitShare)),  // local shares of bmins
            *rbmins= (BitShare  *) malloc(batch_size*sizeof(BitShare)); // remote shares of bmins
  assert(bmins!=NULL); assert(rbmins!=NULL);
  BitShare  *bmaxs= (BitShare  *) malloc(batch_size*sizeof(BitShare)),  // local shares of bmaxs
            *rbmaxs= (BitShare  *) malloc(batch_size*sizeof(BitShare)); // remote shares of bmaxs
  assert(bmaxs!=NULL); assert(rbmaxs!=NULL);
  int len = table->numCols/2 + 2;  // Number of attributes plus new min plus new max
  BShare  *local= (BShare  *) malloc(batch_size*len*sizeof(BShare)),  // local shares of multiplexed rows
          *remote= (BShare  *) malloc(batch_size*len*sizeof(BShare)); // remote shares of multiplexed rows
  assert(local!=NULL); assert(remote!=NULL);
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int rounds = log2(table->numRows);  // Number of rounds for odd-even aggregation
  int dist = table->numRows;  // Distance between elements to compare at each round
  int num_pairs, start, step, end, offset, total_comparisons, num_comparisons=0;
  // Use an odd-even aggregation and count rows by adding the 'selected' bits
  for (int r=0; r<rounds; r++){
    // if (get_rank()==0)
    // printf("Round %d\n", r);

    // Data t[32], ot[32];
    // for (int u=0;u<8; u++) {
    //   t[u*4] = c[u][0];
    //   t[u*4+1] = c[u][2];
    //   t[u*4+2] = c[u][4];
    //   t[u*4+3] = c[u][6];
    // }
    // open_b_array(t, 32, ot);
    // if (get_rank()==0) {
    //   printf("Table at round %d: \n", r);
    //   for (int u=0;u<32;u+=4) {
    //     printf("%lld %lld %lld %lld\n", ot[u], ot[u+1], ot[u+2], ot[u+3]);
    //   }
    //   printf("\n");
    // }

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

      // Data bss[8];
      // open_b_array(local_bits, num_pairs, bss);
      // if (get_rank()==0) {
      //   printf("BSs of round %d: \n", r);
      //   for (int u=0;u<num_pairs;u++) {
      //     printf("%lld ", bss[u]);
      //   }
      //   printf("\n");
      // }

      // if (get_rank()==0)
      // printf("Computing bu.\n");
      // Compute bu = (bs AND selected_b[i] AND selected_b[i+dist]) ("same group and both selected")
      for (int i=start, k=0; i<end; i++, k++) {
        bits[k] = local_bits[k];
        rbits[k] = remote_bits[k];
        bits[k] = and_b(bits[k], rbits[k], 
                        selected_b[i], remote_selected_b[i],
                        get_next_rb());
      }
      exchange_shares_array(bits, rbits, num_pairs);  // 1 round

      // Data bus[8];
      // open_b_array(bits, num_pairs, bus);
      // if (get_rank()==0) {
      //   printf("BUs of round %d: \n", r);
      //   for (int u=0;u<num_pairs;u++) {
      //     printf("%lld ", bus[u]);
      //   }
      //   printf("\n");
      // }

      for (int i=start, k=0; i<end; i++, k++) {
        bits[k] = and_b(bits[k], rbits[k], 
                        selected_b[i+dist], remote_selected_b[i+dist],
                        get_next_rb()) & mask;
      }
      exchange_shares_array(bits, rbits, num_pairs);  // 1 round
      // if (get_rank()==0)
      // printf("Computing bmin.\n");
      // Compute bmin = c[i+dist][min_att] <? c[i][min_att]
      greater_batch2(c, min_att, start, num_pairs, dist, bmins);
      exchange_bit_shares_array(bmins, rbmins, num_pairs);  // 1 round

      // bool omins[8];
      // open_bit_array(bmins, num_pairs, omins);
      // if (get_rank()==0) {
      //   printf("Mins of round %d: \n", r);
      //   for (int u=0;u<num_pairs;u++) {
      //     printf("%d ", omins[u]);
      //   }
      //   printf("\n");
      // }

      // if (get_rank()==0)
      // printf("Computing bmax.\n");
      // Compute bm = c[i+dist][max_att] <? c[i][max_att]
      greater_batch2(c, max_att, start, num_pairs, dist, bmaxs);
      for (int k=0; k<num_pairs; k++) {
        bmaxs[k] ^= mask;  // bmax = NOT bm
      }
      exchange_bit_shares_array(bmaxs, rbmaxs, num_pairs);  // 1 round

      // open_bit_array(bmaxs, num_pairs, omins);
      // if (get_rank()==0) {
      //   printf("Maxs of round %d: \n", r);
      //   for (int u=0;u<num_pairs;u++) {
      //     printf("%d ", omins[u]);
      //   }
      //   printf("\n");
      // }

      // if (get_rank()==0)
      // printf("Computing composite bits.\n");
      // Compute bits:
      //  - bcn = bu AND bmin 
      //  - bcx = bu AND bmax
      for (int k=0; k<num_pairs; k++) {
        bcns[k] = and_b(bits[k], rbits[k], bmins[k], rbmins[k],
                       get_next_rb()) & mask;
        bits[k] = and_b(bits[k], rbits[k], bmaxs[k], rbmaxs[k],
                        get_next_rb()) & mask;
      }
      exchange_shares_array(bcns, rbcns, num_pairs);  // 1 round

      // Data ocomp[8];
      // open_b_array(bcns, num_pairs, ocomp);
      // if (get_rank()==0) {
      //   printf("bcns of round %d: \n", r);
      //   for (int u=0;u<num_pairs;u++) {
      //     printf("%lld ", ocomp[u]);
      //   }
      //   printf("\n");
      // }

      exchange_shares_array(bits, rbits, num_pairs);  // 1 round

      // open_b_array(bits, num_pairs, ocomp);
      // if (get_rank()==0) {
      //   printf("bcxs of round %d: \n", r);
      //   for (int u=0;u<num_pairs;u++) {
      //     printf("%lld ", ocomp[u]);
      //   }
      //   printf("\n");
      // }

      // if (get_rank()==0)
      // printf("Computing bn.\n");
      // Compute masking bit bn = (bs OR NOT selected_b[i+dist]) = NOT(NOT bs AND selected_b[i+dist])
      for (int i=start, k=0; i<end; i++, k++) {
        local_bits[k] ^= mask;
        remote_bits[k] ^= mask; 
        local_bits[k] = and_b(local_bits[k], remote_bits[k], 
                              selected_b[i+dist], remote_selected_b[i+dist],
                              get_next_rb()) ^ mask;
      }
      exchange_shares_array(local_bits, remote_bits, num_pairs);  // 1 round
      // Compute new_c[i+dist] = bn * dummy_row + (1-bn) * c[i+dist] to multiplex rows
      for (int i=start, k=0; i<end; i++, k++) {
        local_bits[k] &= mask;  // Keep LSB only
        local_bits[k] = -local_bits[k];  // Set all bits equal to the LSB
        remote_bits[k] &= mask;  // Keep LSB only
        remote_bits[k] = -remote_bits[k];  // Set all bits equal to the LSB
        offset = k*len;
        for (int j=0; j<table->numCols; j+=2) {
          local[offset+j/2] = and_b(local_bits[k], remote_bits[k], 
                                    max, max, get_next_rb());
          local[offset+j/2] ^= and_b(~local_bits[k], ~remote_bits[k], 
                                     c[i+dist][j], c[i+dist][j+1], get_next_rb());
        }
      }
      for (int i=start, k=0; i<end; i++, k++) {
        offset = (k+1) * len;
        bcns[k] = -bcns[k];
        rbcns[k] = -rbcns[k];
        bits[k] = -bits[k];
        rbits[k] = -rbits[k];
        // Compute new_min = bcn*c[i+dist] + (1-bcn)*c[i]
        local[offset-2] = and_b(bcns[k], rbcns[k], c[i+dist][min_att], c[i+dist][min_att+1],
                                get_next_rb());
        local[offset-2] ^= and_b(~bcns[k], ~rbcns[k], c[i][min_att],
                                 c[i][min_att+1], get_next_rb());
        // Compute new_max = bcx*c[i+dist] + (1-bcx)*c[i]
        local[offset-1] = and_b(bits[k], rbits[k], c[i+dist][max_att], c[i+dist][max_att+1],
                                get_next_rb());
        local[offset-1] ^= and_b(~bits[k], ~rbits[k], c[i][max_att],
                                 c[i][max_att+1], get_next_rb());
      }
      // Fetch remote boolean and arithmetic shares in a single round
      // NOTE: This works because BShare and AShare are both of the same type
      exchange_shares_array(local, remote, num_pairs*len);    // 1 round
      // Update rows in place
      for (int i=start, k=0; i<end; i++, k++) {
        offset = k*len;
        // Set c[i+dist] = new_c[i+dist]
        for (int j=0; j<table->numCols; j+=2) {
          c[i+dist][j] = local[offset+j/2];
          c[i+dist][j+1] = remote[offset+j/2];
        }
        // Set c[i][min_att] = new_min
        c[i][min_att] = local[offset+len-2];
        c[i][min_att+1] = remote[offset+len-2];
        // Set c[i][max_att] = new_max
        c[i][max_att] = local[offset+len-1];
        c[i][max_att+1] = remote[offset+len-1];
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
  int rindex = num_comparisons;  // The index of the first available random number
  // Compute bit = selected[0]==0 to multiplex first row
  BShare bit, rbit;
  local_bits[0] = selected_b[0] ^ 0 ^ (~(BShare)0);
  remote_bits[0] = remote_selected_b[0] ^ 0 ^ (~(BShare)0);
  for (int start=0; start<table->numRows; start+=batch_size) {
    int end = (start+batch_size <= table->numRows ? start+batch_size : table->numRows);
    int comp = end - start;
    offset = start==0 ? 1 : 0;
    // Compute b = c[start][0]==max to multiplex rows and counts
    for (int i=start+offset, k=offset; i<end; i++, k++) {
      local_bits[k] = c[i][0] ^ max ^ (~(BShare)0);
      remote_bits[k] = c[i][1] ^ max ^ (~(BShare)0);
    }
    // Compute equality bits in bulk
    eq_bulk(numlevels, numbits, local_bits, remote_bits, comp); // numlevels-1 rounds
    exchange_shares_array(local_bits, remote_bits, comp);  // 1 round

    for (int i=start, k=0; i<end; i++, k++) {
      bit = local_bits[k] & mask;  // Keep LSB only
      bit = -bit;  // Set all bits equal to the LSB
      rbit = remote_bits[k] & mask;  // Keep LSB only
      rbit = -rbit;  // Set all bits equal to the LSB
      offset = k*len;
      // Compute c[i] = bit*max + (1-bit)*c[i]
      for (int j=0; j<table->numCols; j+=2) {
        local[offset+j/2] = and_b(bit, rbit, max, max, get_next_rb());
        local[offset+j/2] ^= and_b(~bit, ~rbit, c[i][j], c[i][j+1], get_next_rb());
      }
      offset += len;
      // Compute selected[i] = b*max + (1-b)*selected_b[i]
      local[offset-2] = and_b(bit, rbit, max, max, get_next_rb());
      local[offset-2] ^= and_b(~bit, ~rbit, selected_b[i], remote_selected_b[i],
                               get_next_rb());
    }
    exchange_shares_array(local, remote, comp*len); // 1 round
    rindex += comp;
    // Update rows and counts
    for (int i=start, k=0; i<end; i++, k++) {
      offset = k*len;
      for (int j=0; j<table->numCols; j+=2) {
        c[i][j] = local[offset+j/2];
        c[i][j+1] = remote[offset+j/2];
      }
      offset += len;
      selected_b[i] = local[offset-2];
      remote_selected_b[i] = remote[offset-2];
    }
  }
  // NOTE: Final shuffle iff not followed by ORDER_BY
  free(remote_selected_b);
  free(bits); free(rbits);
  free(bcns); free(rbcns);
  free(bmins); free(rbmins);
  free(bmaxs); free(rbmaxs);
  free(local_bits); free(remote_bits);
  free(local); free(remote);
}

// Groups rows and counts the number of rows per group
// The given BShareTable must be sorted on the group-by key(s)
void group_by_count_micro(BShareTable* table, unsigned* key_indices,
                          int num_keys, AShare* counters,
                          AShare* remote_counters, BShare* rb, AShare* ra) {
  BShare** c = table->content;
  BShare max=0xFFFFFFFFFFFFFFFF;

  // Scan table and update counts by adding 'selected' bits
  AShare ab1, ab2;
  int len = table->numCols/2 + 2;   // Number of attributes plus new and old cnt
  int rank = get_rank();
  int succ_rank = get_succ();
  BShare bs1, bs2;
  BShare local[len], remote[len];
  BShare local_bits[num_keys], remote_bits[num_keys];
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
    BShare b1 = -bs1;                   // Set all bits equal to LSB of bs1
    BShare b2 = -bs2;                   // Set all bits equal to LSB of bs2
    // Compute new_c[i] = b * dummy_row + (1-b) * row1
    for (int j=0; j<table->numCols; j+=2) {
      local[j/2] = and_b(b1, b2, max, max, get_next_rb());
      local[j/2] ^= and_b(~b1, ~b2, c[i][j], c[i][j+1], get_next_rb());
    }
    // Compute arithmetic shares from boolean shares
    ab1 = convert_single_bit(bs1, ra[k], rb[k]);    // 1 round
    ab2 = exchange_shares(ab1);                     // 1 round
    // Compute new_cnt = bs*(counters[i] + counters[i+1]) + (1-bs)*counters[i+1]
    AShare local_count = counters[i] + counters[i+1],
           remote_count = remote_counters[i] + remote_counters[i+1];
    local[len-2] = mul(ab1, ab2, local_count, remote_count,
                       get_next_r());
    local[len-2] += mul(rank%2 - ab1, succ_rank%2 - ab2, counters[i+1],
                        remote_counters[i+1], get_next_r());
    // Compute ab * max + (1-ab) * counters[i] to allow masking previous count
    local[len-1] = mul(ab1, ab2, max, max, get_next_r());
    local[len-1] += mul(rank%2 - ab1, succ_rank%2 - ab2,
                        counters[i], remote_counters[i], get_next_r());
    // Fetch remote boolean and arithmetic shares
    // NOTE: This works because BShare and AShare are both long long ints
    exchange_shares_array(local, remote, len);    // 1 round
    // Set c[i] = new_c[i]
    for (int j=0; j<table->numCols; j+=2) {
      c[i][j] = local[j/2];
      c[i][j+1] = remote[j/2];
    }
    // Update 'selected' bit
    counters[i] = local[len-1];
    remote_counters[i] = remote[len-1];
    // Set selected[i+1] = new_cnt
    counters[i+1] = local[len-2];
    remote_counters[i+1] = remote[len-2];
  }
  // NOTE: Final shuffle iff not followed by ORDER_BY
}

// Fused group-by-join-aggregattion (sum) operation
// on two input tables, left and right.
void group_by_join(BShareTable* left, BShareTable* right, int start_left,
                   int end_left, int group_att_index, int left_join_index,
                   int right_join_index, int right_att_index, BShare* rb_left,
                   AShare* ra_left, BShare* rb_right, AShare* ra_right,
                   unsigned sum_res_index, unsigned count_res_index) {

  int rank = get_rank();
  int succ_rank = get_succ();

  int batch_size = end_left - start_left;

  // Equality bits
  BShare *b2 = (BShare *) malloc(batch_size*right->numRows*sizeof(BShare));
  assert(b2!=NULL);
  // Remote equality bits
  BShare *b2_remote = (BShare *) malloc(batch_size*right->numRows*sizeof(BShare));
  assert(b2_remote!=NULL);
  // Converted b2
  AShare *b2_a = (AShare *) malloc(batch_size*right->numRows*sizeof(AShare));
  assert(b2_a!=NULL);
  // Partial aggregate
  AShare *sum_right = (AShare *) calloc(batch_size, sizeof(AShare));
  assert(sum_right!=NULL);

  // Step 1: compute the equality predicate for one row of the left input
  // and all rows on the right input in a batch
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int index;
  // For each row in the batch
  for (int i=start_left, p=0; i<end_left; i++, p++) {
    // Compute bitwise x^y^1
    for (int j=0; j<right->numRows; j++) {
      index = p*right->numRows + j;
      b2[index] = left->content[i][left_join_index] ^
                  right->content[j][right_join_index] ^
                  (~(BShare)0); // local share;
      b2_remote[index] = left->content[i][left_join_index + 1] ^
                         right->content[j][right_join_index + 1] ^
                         (~(BShare)0); // remote share
    }
  }
  // The result is stored in the (numbits/2) rightmost bits of result
  for (int l=0; l<numlevels; l++) {
    // For each row in the batch
    for (int i=start_left, p=0; i<end_left; i++, p++) {
      // For each row on the right
      for (int k=0; k<right->numRows; k++) {
        index = p*right->numRows + k;
        b2[index] = eq_b_level2(numbits >> l, b2[index], b2_remote[index]);
      }
    }
    // Exchange results of logical and, except for the final round
    if (l != numlevels-1) {
      exchange_shares_array(b2, b2_remote, batch_size*right->numRows);
    }
  }

  // Step 2: Convert equality bits to arithmetic shares
  convert_single_bit_array(b2, ra_right, rb_right,
                            batch_size*right->numRows, b2_a);
  // exchange arithmetic bits
  exchange_shares_array(b2_a, b2_remote, batch_size*right->numRows);

  // Step 3: Compute the row sum as s+=b2[j] * score[j]
  // For each row in the batch
  for (int i=start_left, p=0; i<end_left; i++, p++) {
    for (int k=0; k<right->numRows; k++) {
      index = p*right->numRows + k;
      sum_right[p] += mul(b2_a[index], b2_remote[index],
                      right->content[k][right_att_index],
                      right->content[k][right_att_index + 1],
                      get_next_r());
    }
  }
  // Free memory
  free(b2); free(b2_a); free(b2_remote);

  // Remote partial aggregate
  AShare *sum_right_remote = (AShare *) calloc(batch_size, sizeof(AShare));
  assert(sum_right_remote!=NULL);

  // Get remote partial sums
  exchange_shares_array(sum_right, sum_right_remote, batch_size);
  // Group-by bits
  BShare *b1 = (BShare *) malloc(batch_size*sizeof(BShare));
  assert(b1!=NULL);
  // Group-by remote bits
  BShare *b1_remote = (BShare *) malloc(batch_size*sizeof(BShare));
  assert(b1_remote!=NULL);
  // Converted group-by bits
  AShare *b1_a = (AShare *) malloc(batch_size*sizeof(AShare));
  assert(b1_a!=NULL);
  // Converted remote group-by bits
  AShare *b1_a_remote = (AShare *) malloc(batch_size*sizeof(AShare));
  assert(b1_a_remote!=NULL);

  // Step 4: Compute b1 (group-by), convert it to arithmetic and then
  distinct_batch_incr(left, start_left, end_left, group_att_index, b1);
  // Exchange boolean bits
  exchange_shares_array(b1, b1_remote, batch_size);
  // Convert equality bits to arithmetic shares
  convert_single_bit_array(b1, ra_left, rb_left, batch_size, b1_a);
  // Exchange arithmetic bits
  exchange_shares_array(b1_a, b1_a_remote, batch_size);

  // Aggregate and mask
  BShare max=0xFFFFFFFFFFFFFFFF;
  index = start_left;
  int pos = 0;
  if (index==0) { // If it's the very first row
    left->content[index][sum_res_index] = sum_right[0];
    left->content[index][sum_res_index + 1] = sum_right_remote[0];
    index += 1; // Start aggregating from the second row
    pos += 1;
  }

  for (int i=index; i<end_left; i++) {
    // Compute res = (1-b1)*(prev + sum_right) + b1*sum_right
    left->content[i][sum_res_index] = mul(rank % 2 - b1_a[pos],
                     succ_rank%2 - b1_a_remote[pos],
                                          left->content[i - 1][sum_res_index] + sum_right[pos],
                                          left->content[i - 1][sum_res_index + 1] + sum_right_remote[pos],
                                          get_next_r());

    left->content[i][sum_res_index] += mul(b1_a[pos], b1_a_remote[pos],
                                           sum_right[pos], sum_right_remote[pos],
                                           get_next_r());
    // Exchange the result aggregation to use it in the next iteration
    exchange_shares_array(&left->content[i][sum_res_index],
                          &left->content[i][sum_res_index + 1], 1);
    BShare bb1 = -b1[pos];
    BShare bb2 = -b1_remote[pos];
    // Compute row_{i-1} = (1-b1)*max + b1*row_{i-1}
    for (int j=0; j<left->numCols-1; j+=2) {
      // NOTE: Some attributes are boolean shares some others arithmetic
      if ( (j==sum_res_index) ||
           (j==count_res_index) ) { // It's an arithmetic share
        AShare left_att = left->content[i - 1][j];
        AShare left_att_rem = left->content[i - 1][j + 1];
        left->content[i - 1][j] = mul(rank % 2 - b1_a[pos],
                                     succ_rank%2 - b1_a_remote[pos],
                                      max,
                                      max,
                                      get_next_r());
        left->content[i - 1][j] += mul(b1_a[pos], b1_a_remote[pos],
                                       left_att,
                                       left_att_rem,
                                       get_next_r());
      }
      else {  // It's a boolean share
        BShare left_att = left->content[i - 1][j];
        BShare left_att_rem = left->content[i - 1][j + 1];
        left->content[i - 1][j] = and_b(~bb1, ~bb2,
                                        max,
                                        max,
                                        get_next_rb());
        left->content[i - 1][j] ^= and_b(bb1, bb2,
                                         left_att,
                                         left_att_rem,
                                         get_next_rb());
      }
    }
    pos++;
  }
  // Free memory
  free(b1); free(b1_a); free(b1_remote); free(b1_a_remote);
  free(sum_right); free(sum_right_remote);
}

// Applies first phase of join-aggregation decomposition
void group_by_join_first(BShareTable* left, BShareTable* right, int start_left,
                            int end_left, int left_join_index,
                            int right_join_index, int right_att_index,
                            BShare* rb_right, AShare* ra_right,
                            unsigned sum_res_index) {

  // Make sure the end is at most numRows
  end_left = end_left > left->numRows ? left->numRows : end_left;

  int batch_size = end_left - start_left;

  // Equality bits
  BShare *b2 = (BShare *) malloc(batch_size*right->numRows*sizeof(BShare));
  assert(b2!=NULL);
  // Remote equality bits
  BShare *b2_remote = (BShare *) malloc(batch_size*right->numRows*sizeof(BShare));
  assert(b2_remote!=NULL);
  // Converted b2
  AShare *b2_a = (AShare *) malloc(batch_size*right->numRows*sizeof(AShare));
  assert(b2_a!=NULL);
  // Partial aggregate
  AShare *sum_right = (AShare *) calloc(batch_size, sizeof(AShare));
  assert(sum_right!=NULL);
  // Remote partial aggregate
  AShare *sum_right_remote = (AShare *) calloc(batch_size, sizeof(AShare));
  assert(sum_right_remote!=NULL);
  // Step 1: compute the equality predicate for one row of the left input
  // and all rows on the right input in a batch
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int index;
  // For each row in the batch
  for (int i=start_left, p=0; i<end_left; i++, p++) {
    // Compute bitwise x^y^1
    for (int j=0; j<right->numRows; j++) {
      index = p*right->numRows + j;
      b2[index] = left->content[i][left_join_index] ^
                  right->content[j][right_join_index] ^
                  (~(BShare)0); // local share;
      b2_remote[index] = left->content[i][left_join_index + 1] ^
                         right->content[j][right_join_index + 1] ^
                         (~(BShare)0); // remote share
    }
  }
  // The result is stored in the (numbits/2) rightmost bits of result
  for (int l=0; l<numlevels; l++) {
    // For each row in the batch
    for (int i=start_left, p=0; i<end_left; i++, p++) {
      // For each row on the right
      for (int k=0; k<right->numRows; k++) {
        index = p*right->numRows + k;
        b2[index] = eq_b_level2(numbits >> l, b2[index], b2_remote[index]);
      }
    }
    // Exchange results of logical and, except for the final round
    if (l != numlevels-1) {
      exchange_shares_array(b2, b2_remote, batch_size*right->numRows);
    }
  }

  // Step 2: Convert equality bits to arithmetic shares
  convert_single_bit_array(b2, ra_right, rb_right,
                           batch_size*right->numRows, b2_a);
  // exchange arithmetic bits
  exchange_shares_array(b2_a, b2_remote, batch_size*right->numRows);

  // Step 3: Compute the row sum as s+=b2[j] * score[j]
  for (int i=start_left, p=0; i<end_left; i++, p++) {
    for (int k=0; k<right->numRows; k++) {
      index = p*right->numRows + k;
      sum_right[p] += mul(b2_a[index], b2_remote[index],
                      right->content[k][right_att_index],
                      right->content[k][right_att_index + 1],
                      get_next_r());
    }
  }
  // Get remote partial sums
  exchange_shares_array(sum_right, sum_right_remote, batch_size);
  // Set partial sums
  for (int i=start_left, p=0; i<end_left; i++, p++) {
    left->content[i][sum_res_index] = sum_right[p];
    left->content[i][sum_res_index + 1] = sum_right_remote[p];
  }
  // Free memory
  free(b2); free(b2_a); free(b2_remote);
  free(sum_right); free(sum_right_remote);
}

// Applies second phase of join-aggregation decomposition
void group_by_sum_odd_even(BShareTable* left, int batch_size,
                           BShare* rb_left, AShare* ra_left, unsigned sum_res_index,
                           unsigned count_res_index, unsigned* key_indices, unsigned num_keys) {
  // Make sure the number of rows is a power of two
  assert(ceil(log2(left->numRows)) == floor(log2(left->numRows)));
  // Make sure batch_size is at most equal to the input size
  batch_size = batch_size > left->numRows ? left->numRows : batch_size;
  BShare** c = left->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Allocate memory for local and remote shares
  BShare  *local_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)),  // local shares of equality bits
          *remote_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)); // remote shares of equality bits
  assert(local_bits!=NULL); assert(remote_bits!=NULL);
  AShare  *ab1= (AShare  *) malloc(batch_size*sizeof(AShare)),  // local shares of arithmetic equality bits
          *ab2= (AShare  *) malloc(batch_size*sizeof(AShare));  // remote shares of arithmetic equality bits
  assert(ab1!=NULL); assert(ab2!=NULL);
  int len = left->numCols/2 + 1;  // Number of attributes plus new count
  BShare  *local= (BShare  *) malloc(batch_size*len*sizeof(BShare)),  // local shares of multiplexed rows
          *remote= (BShare  *) malloc(batch_size*len*sizeof(BShare)); // remote shares of multiplexed rows
  assert(local!=NULL); assert(remote!=NULL);
  int rank = get_rank();
  int succ_rank = get_succ();
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int rounds = log2(left->numRows);  // Number of rounds for odd-even aggregation
  int dist = left->numRows;  // Distance between elements to compare at each round
  int num_pairs, start, step, end, offset, total_comparisons, num_comparisons=0;
  // Use an odd-even aggregation and count rows by adding the 'selected' bits
  for (int r=0; r<rounds; r++){
    dist /= 2;
    total_comparisons = left->numRows-dist;
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
      // Compute arithmetic shares of the equality bits using the boolean shares
      convert_single_bit_array(local_bits, &ra_left[num_comparisons], &rb_left[num_comparisons], 
                               num_pairs, ab1);   // 2 rounds
      exchange_shares_array(ab1, ab2, num_pairs); // 1 round
      // Compute new_c[i+dist] = bs*dummy_row + (1-bs)*c[i+dist] to multiplex rows
      for (int i=start, k=0; i<end; i++, k++) {
        local_bits[k] &= mask;  // Keep LSB only
        local_bits[k] = -local_bits[k];  // Set all bits equal to the LSB
        remote_bits[k] &= mask;  // Keep LSB only
        remote_bits[k] = -remote_bits[k];  // Set all bits equal to the LSB
        offset = k*len;
        AShare a_bs1 = ab1[k], a_bs2 = ab2[k];
        for (int j=0; j<left->numCols; j+=2) {
          if (j==sum_res_index || j==count_res_index) {
            local[offset+j/2] = mul(a_bs1, a_bs2, 0, 0, get_next_r());
            local[offset+j/2] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, c[i+dist][j], c[i+dist][j+1], get_next_r());
          }
          else{
            local[offset+j/2] = and_b(local_bits[k], remote_bits[k], max, max, get_next_rb());
            local[offset+j/2] ^= and_b(~local_bits[k], ~remote_bits[k], c[i+dist][j], c[i+dist][j+1], get_next_rb());
          }
        }
        // Compute new_cnt = bs*(counters[i] + counters[i+dist]) + (1-bs)*counters[i]
        AShare local_count = c[i][sum_res_index] + c[i+dist][sum_res_index],
               remote_count = c[i][sum_res_index+1] + c[i+dist][sum_res_index+1];
        offset = (k+1) * len;
        local[offset-1] = mul(a_bs1, a_bs2, local_count, remote_count,
                              get_next_r());
        local[offset-1] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, c[i][sum_res_index],
                               c[i][sum_res_index+1], get_next_r());
      }
      // Fetch remote boolean and arithmetic shares in a single round
      // NOTE: This works because BShare and AShare are both of the same type
      exchange_shares_array(local, remote, num_pairs*len);    // 1 round
      // Update rows in place
      for (int i=start, k=0; i<end; i++, k++) {
        offset = k*len;
        // Set c[i+dist] = new_c[i+dist]
        for (int j=0; j<left->numCols; j+=2) {
          c[i+dist][j] = local[offset+j/2];
          c[i+dist][j+1] = remote[offset+j/2];
        }
        // Set c[i][sum_res_index] = new_cnt
        c[i][sum_res_index] = local[offset+len-1];
        c[i][sum_res_index+1] = remote[offset+len-1];
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
  int rindex = num_comparisons;
  for (int start=0; start<left->numRows; start+=batch_size) {
    int end = (start+batch_size <= left->numRows ? start+batch_size : left->numRows);
    int comp = end - start;
    // Compute b = c[start][0]==max that denotes whether the row has been masked (1) or not (0)
    for (int i=start, k=0; i<end; i++, k++) {
      local_bits[k] = c[i][0] ^ max ^ (~(BShare)0); 
      remote_bits[k] = c[i][1] ^ max ^ (~(BShare)0);
    }
    // Compute equality bits in bulk
    eq_bulk(numlevels, numbits, local_bits, remote_bits, comp); // numlevels-1 rounds
    exchange_shares_array(local_bits, remote_bits, comp);  // 1 round
    // Compute arithmetic shares of the equality bits using the boolean shares
    convert_single_bit_array(local_bits, &ra_left[rindex], &rb_left[rindex], 
                             comp, ab1);   // 2 rounds
    exchange_shares_array(ab1, ab2, comp); // 1 round
    rindex += comp;
    // Compute c[i] = b*max + (1-b)*c[i]
    for (int i=start, k=0; i<end; i++, k++) {
      local_bits[k] &= mask;
      remote_bits[k] &= mask;
      local_bits[k] = -local_bits[k];
      remote_bits[k] = -remote_bits[k];
      offset = k * len;
      AShare a_bs1 = ab1[k], a_bs2 = ab2[k];
      for (int j=0; j<left->numCols; j+=2) {
        if (j==sum_res_index || j==count_res_index) {
            local[offset+j/2] = mul(a_bs1, a_bs2, 0, 0, get_next_r());
            local[offset+j/2] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, c[i][j], c[i][j+1], get_next_r());
        }
        else {
          local[offset+j/2] = and_b(local_bits[k], remote_bits[k], max, max, get_next_rb());
          local[offset+j/2] ^= and_b(~local_bits[k], ~remote_bits[k],
                                     c[i][j], c[i][j+1], get_next_rb());
        }  
      }
    }
    exchange_shares_array(local, remote, comp*len); // 1 round
    // Update rows
    for (int i=start, k=0; i<end; i++, k++) {
      offset = k * len;
      for (int j=0; j<left->numCols; j+=2) {
        c[i][j] = local[offset+j/2];
        c[i][j+1] = remote[offset+j/2];
      }      
    }
  }
  // NOTE: Final shuffle iff not followed by ORDER_BY
  free(local_bits); free(remote_bits);
  free(ab1); free(ab2);
  free(local); free(remote); 
}

// Computes the average (sum, count) and updates relation in place
// Expects relation to be sorted on the group-by keys (key_indices)
// Expects column at count_res_index to be initiailized to `1`
void group_by_avg_odd_even(BShareTable* rel, int batch_size,
                           BShare* rb_rel, AShare* ra_rel, unsigned sum_res_index,
                           unsigned count_res_index, unsigned* key_indices, unsigned num_keys) {
  // Make sure the number of rows is a power of two
  assert(ceil(log2(rel->numRows)) == floor(log2(rel->numRows)));
  // Make sure batch_size is at most equal to the input size
  batch_size = batch_size > rel->numRows ? rel->numRows : batch_size;
  BShare** c = rel->content;
  BShare mask=1, max=0xFFFFFFFFFFFFFFFF;
  // Allocate memory for local and remote shares
  BShare  *local_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)),  // local shares of equality bits
          *remote_bits= (BShare  *) malloc(batch_size*num_keys*sizeof(BShare)); // remote shares of equality bits
  assert(local_bits!=NULL); assert(remote_bits!=NULL);
  AShare  *ab1= (AShare  *) malloc(batch_size*sizeof(AShare)),  // local shares of arithmetic equality bits
          *ab2= (AShare  *) malloc(batch_size*sizeof(AShare));  // remote shares of arithmetic equality bits 
  assert(ab1!=NULL); assert(ab2!=NULL);
  int len = rel->numCols/2 + 2;  // Number of attributes plus new sum plus new count
  BShare  *local= (BShare *) malloc(batch_size*len*sizeof(BShare)),  // local shares of multiplexed rows
          *remote= (BShare *) malloc(batch_size*len*sizeof(BShare)); // remote shares of multiplexed rows
  assert(local!=NULL); assert(remote!=NULL);
  int rank = get_rank();
  int succ_rank = get_succ();
  int numbits = sizeof(BShare) * 8;
  int numlevels = log2(numbits);
  int rounds = log2(rel->numRows);  // Number of rounds for odd-even aggregation
  int dist = rel->numRows;  // Distance between elements to compare at each round
  int num_pairs, start, step, end, offset, total_comparisons, num_comparisons=0;
  // Use an odd-even aggregation and count rows by adding the 'selected' bits
  for (int r=0; r<rounds; r++){
    dist /= 2;
    total_comparisons = rel->numRows-dist;
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
      // Compute arithmetic shares of the equality bits using the boolean shares
      convert_single_bit_array(local_bits, &ra_rel[num_comparisons], &rb_rel[num_comparisons], 
                               num_pairs, ab1);   // 2 rounds
      exchange_shares_array(ab1, ab2, num_pairs); // 1 round
      // Compute new_c[i+dist] = bs*dummy_row + (1-bs)*c[i+dist] to multiplex rows
      for (int i=start, k=0; i<end; i++, k++) {
        local_bits[k] &= mask;  // Keep LSB only
        local_bits[k] = -local_bits[k];  // Set all bits equal to the LSB
        remote_bits[k] &= mask;  // Keep LSB only
        remote_bits[k] = -remote_bits[k];  // Set all bits equal to the LSB
        offset = k*len;
        AShare a_bs1 = ab1[k], a_bs2 = ab2[k];
        for (int j=0; j<rel->numCols; j+=2) {
          if (j==sum_res_index || j==count_res_index) {
            local[offset+j/2] = mul(a_bs1, a_bs2, 0, 0, get_next_r());
            local[offset+j/2] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, c[i+dist][j], c[i+dist][j+1], get_next_r());
          }
          else{
            local[offset+j/2] = and_b(local_bits[k], remote_bits[k], max, max, get_next_rb());
            local[offset+j/2] ^= and_b(~local_bits[k], ~remote_bits[k], c[i+dist][j], c[i+dist][j+1], get_next_rb());
          }
        }
        // Compute new_sum = bs*(c[i] + c[i+dist]) + (1-bs)*c[i]
        AShare local_sum = c[i][sum_res_index] + c[i+dist][sum_res_index],
               remote_sum = c[i][sum_res_index+1] + c[i+dist][sum_res_index+1];
        offset = (k+1) * len;
        local[offset-2] = mul(a_bs1, a_bs2, local_sum, remote_sum,
                              get_next_r());
        local[offset-2] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, c[i][sum_res_index],
                               c[i][sum_res_index+1], get_next_r());
        // Compute new_cnt = bs*(c[i] + c[i+dist]) + (1-bs)*c[i]
        AShare local_count = c[i][count_res_index] + c[i+dist][count_res_index],
               remote_count = c[i][count_res_index+1] + c[i+dist][count_res_index+1];
        local[offset-1] = mul(a_bs1, a_bs2, local_count, remote_count,
                              get_next_r());
        local[offset-1] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, c[i][count_res_index],
                               c[i][count_res_index+1], get_next_r());
      }
      // Fetch remote boolean and arithmetic shares in a single round
      // NOTE: This works because BShare and AShare are both of the same type
      exchange_shares_array(local, remote, num_pairs*len);    // 1 round
      // Update rows in place
      for (int i=start, k=0; i<end; i++, k++) {
        offset = k*len;
        // Set c[i+dist] = new_c[i+dist]
        for (int j=0; j<rel->numCols; j+=2) {
          c[i+dist][j] = local[offset+j/2];
          c[i+dist][j+1] = remote[offset+j/2];
        }
        // Set c[i][sum_res_index] = new_sum
        c[i][sum_res_index] = local[offset+len-2];
        c[i][sum_res_index+1] = remote[offset+len-2];
        // Set c[i][sum_count_index] = new_cnt
        c[i][count_res_index] = local[offset+len-1];
        c[i][count_res_index+1] = remote[offset+len-1];
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
  int rindex = num_comparisons;
  for (int start=0; start<rel->numRows; start+=batch_size) {
    int end = (start+batch_size <= rel->numRows ? start+batch_size : rel->numRows);
    int comp = end - start;
    // Compute b = c[start][0]==max that denotes whether the row has been masked (1) or not (0)
    for (int i=start, k=0; i<end; i++, k++) {
      local_bits[k] = c[i][0] ^ max ^ (~(BShare)0); 
      remote_bits[k] = c[i][1] ^ max ^ (~(BShare)0);
    }
    // Compute equality bits in bulk
    eq_bulk(numlevels, numbits, local_bits, remote_bits, comp); // numlevels-1 rounds
    exchange_shares_array(local_bits, remote_bits, comp);  // 1 round
    // Compute arithmetic shares of the equality bits using the boolean shares
    convert_single_bit_array(local_bits, &ra_rel[rindex], &rb_rel[rindex], 
                             comp, ab1);   // 2 rounds
    exchange_shares_array(ab1, ab2, comp); // 1 round
    rindex += comp;
    // Compute c[i] = b*max + (1-b)*c[i]
    for (int i=start, k=0; i<end; i++, k++) {
      local_bits[k] &= mask;
      remote_bits[k] &= mask;
      local_bits[k] = -local_bits[k];
      remote_bits[k] = -remote_bits[k];
      offset = k * len;
      AShare a_bs1 = ab1[k], a_bs2 = ab2[k];
      for (int j=0; j<rel->numCols; j+=2) {
        if (j==sum_res_index || j==count_res_index) {
            local[offset+j/2] = mul(a_bs1, a_bs2, 0, 0, get_next_r());
            local[offset+j/2] += mul(rank%2 - a_bs1, succ_rank%2 - a_bs2, c[i][j], c[i][j+1], get_next_r());
        }
        else {
          local[offset+j/2] = and_b(local_bits[k], remote_bits[k], max, max, get_next_rb());
          local[offset+j/2] ^= and_b(~local_bits[k], ~remote_bits[k],
                                     c[i][j], c[i][j+1], get_next_rb());
        }  
      }
    }
    exchange_shares_array(local, remote, comp*len); // 1 round
    // Update rows
    for (int i=start, k=0; i<end; i++, k++) {
      offset = k * len;
      for (int j=0; j<rel->numCols; j+=2) {
        c[i][j] = local[offset+j/2];
        c[i][j+1] = remote[offset+j/2];
      }      
    }
  }
  // NOTE: Final shuffle iff not followed by ORDER_BY
  free(local_bits); free(remote_bits);
  free(ab1); free(ab2);
  free(local); free(remote); 
}

PRIVATE void bitonic_merge(BShare** contents, int low, int cnt,
                           unsigned index_1, unsigned index_2,
                           int num_elements, int asc) {
  if (cnt>1) {
    int k = cnt/2;
    for (int i=low; i<low+k; i++) {
      // Compare rows i, i+k and swap if necessary
      cmp_swap_g(contents[i], contents[i+k], index_1, index_2, num_elements,
                 asc);
    }
    bitonic_merge(contents, low, k, index_1, index_2, num_elements, asc);
    bitonic_merge(contents, low+k, k, index_1, index_2, num_elements, asc);
  }
}

// Sorts the given table of BShares in place
void bitonic_sort(BShareTable* table, int low, int cnt, unsigned sort_attribute,
                  int asc) {
  if (cnt>1) {
    int k = cnt/2;
    // sort in ascending order since asc here is 1
    bitonic_sort(table, low, k, sort_attribute, 1);
    // sort in descending order since asc here is 0
    bitonic_sort(table, low+k, k, sort_attribute, 0);
    // Will merge whole sequence in ascending order
    // since asc=1
    bitonic_merge(table->content, low, cnt, sort_attribute, sort_attribute,
                  table->numCols, asc);
  }
}

// Same result as bitonic_sort() but works in batch mode
void bitonic_sort_batch(BShareTable* table, unsigned* sort_attributes,
                        int num_attributes, bool* asc, int batch_size) {
  // Batch size and table size must both be a power of two
  assert(ceil(log2(table->numRows)) == floor(log2(table->numRows)));
  assert(ceil(log2(batch_size)) == floor(log2(batch_size)));
  // Batch size must be less than or equal to n/2
  assert(batch_size <= (table->numRows / 2));
  int rounds = (int) log2(table->numRows);
  int num_batches = (table->numRows / 2) / batch_size;
  for (int i = 0; i < rounds; i++) {
    for(int j = 0; j <= i; j++) {
      int last_index = 0;
      for (int r=0; r<num_batches; r++){
        last_index = cmp_swap_batch(last_index, table->content,
                                      table->numRows, table->numCols,
                                      sort_attributes, num_attributes, i, j,
                                      asc, batch_size);
      }
    }
  }
}

// Masks the non-selected rows in the given table
void mask(BShareTable* table, BShare* selected, int batch_size) {
  BShare max=0xFFFFFFFFFFFFFFFF;
  BShare** c = table->content;
  // Fetch the second arithmetic share of each 'selected' bit -- 1 round
  BShare* remote_selected = (BShare* ) malloc((table->numRows)*sizeof(BShare));
  assert(remote_selected!=NULL);
  exchange_shares_array(selected, remote_selected, table->numRows);
  // Make sure batch size is not_b larger than the input table
  batch_size = (batch_size > table->numRows ? table->numRows : batch_size);
  // Allocate batches for local and remote shares
  int width = table->numCols/2, len = width * batch_size;
  BShare* local = (BShare* ) malloc(len*sizeof(BShare));
  assert(local!=NULL);
  BShare* remote = (BShare* ) malloc(len*sizeof(BShare));
  assert(remote!=NULL);
  int start=0, end=batch_size, step;
  BShare b1, b2;
  // For all rows in the input table
  while (start<table->numRows) {
    // For all rows within the given batch size
    for (int i=start, k=0; i<end; i++, k+=width) {
      b1 = selected[i] & 1;
      b2 = remote_selected[i] & 1;
      // Compute c_dummy = selected[i]*c[i] + (1-selected[i])*max
      for (int j=0; j<table->numCols; j+=2) {
        local[k+j/2] = and_b(-b1, -b2, c[i][j], c[i][j+1], get_next_rb());
        local[k+j/2] ^= and_b(~(-b1), ~(-b2), max, max, get_next_rb());
      }
    }
    // Fetch remote shares - 1 round
    exchange_shares_array(local, remote, len);
    for (int i=start, k=0; i<end; i++, k+=width) {
      // Set new row
      for (int j=0; j<table->numCols; j+=2) {
        c[i][j] = local[k+j/2];
        c[i][j+1] = remote[k+j/2];
      }
    }
    start = end;
    step = end + batch_size;
    end = (step <= table->numRows ? step : table->numRows);
  }
  free(remote_selected); free(local); free(remote);
}