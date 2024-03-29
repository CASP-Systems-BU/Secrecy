#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "comm.h"
#include "mpctypes.h"
#include "party.h"
#include "utils.h"


/*******************************************************************
 * add: Adds two numbers by adding their arithmetic shares x, y
 *      Result: x + y
 *      Does not require party-to-party communication
 * *****************************************************************/
AShare add(AShare x, AShare y);

/*******************************************************************
 * mul: Multiplies two numbers x, y using arithmetic shares
 *      Assumes three computing parties holding the shares:
 *      x = x1 + x2 + x3
 *      y = y1 + y2 + y3
 *
 *      P1 has: x1, x2, y1, y2
 *      P2 has: x2, x3, y2, y3
 *      P3 has: x3, x1, y3, y1
 *
 *      The three parties are positioned on a ring. Each party has generated
 *      a seed and sent it to its successor party. Each party ends up with
 *      two seeds: its own and the seed from its predecessor on the ring.
 *      Each party uses these two seeds to generate two random numbers. Then,
 *      each party uses the random numbers to generate its rnum as follows:
 *      P1: P1_rnum = random(P1_seed) - random(P3_seed)
 *      P2: P2_rnum = random(P2_seed) - random(P1_seed)
 *      P3: P3_rnum = random(P3_seed) - random(P2_seed)
 *      The following holds: P1_rnum + P2_rnum + P3_rnum = 0
 *
 *      P1 calls mul() with input: x1, x2, y1, y2, P1_rnum
 *      P2 calls mul() with input: x2, x3, y2, y3, P2_rnum
 *      P3 calls mul() with input: x3, x1, y3, y1, P3_rnum
 *
 *      Does NOT require Beaver triples or all-to-all communication
 *      BUT requires pseudo-random sharing of zero
 * *****************************************************************/
AShare mul(AShare x1, AShare x2, AShare y1, AShare y2, AShare rnum);

/*******************************************************************
 * eq: Checks equality of two numbers x, y using arithmetic shares
 *      x == y <=> (x - y) * w == 0, where w != 0
 *      Assumes three computing parties holding the shares:
 *      x = x1 + x2 + x3
 *      y = y1 + y2 + y3
 *      w = w1 + w2 + w3, each one generated by P1, P2, and P3 respectively
 *
 *      P1 has: x1, x2, y1, y2, w1, w2
 *      P2 has: x2, x3, y2, y3, w2, w3
 *      P3 has: x3, x1, y3, y1, w3, w1
 *
 *      The three parties are positioned on a ring. Each party has generated
 *      a seed and sent it to its successor party. Each party ends up with
 *      two seeds: its own and the seed from its predecessor on the ring.
 *      Each party uses these two seeds to generate two random numbers. Then,
 *      each party uses the random numbers to generate its rnum as follows:
 *      P1: P1_rnum = random(P1_seed) - random(P3_seed)
 *      P2: P2_rnum = random(P2_seed) - random(P1_seed)
 *      P3: P3_rnum = random(P3_seed) - random(P2_seed)
 *      The following holds: P1_rnum + P2_rnum + P3_rnum = 0
 *
 *      P1 calls eq() with input: x1, x2, y1, y2, w1, w2, P1_rnum
 *      P2 calls eq() with input: x2, x3, y2, y3, w2, w3, P2_rnum
 *      P3 calls eq() with input: x3, x1, y3, y1, w3, w1, P3_rnum
 *
 *      Does NOT require Beaver triples or all-to-all communication
 *      BUT requires pseudo-random sharing of zero
 * *****************************************************************/
AShare eq(AShare x1, AShare x2, AShare y1, AShare y2, AShare w1, AShare w2,
         AShare rnum);

/*******************************************************************
* or_a: Computes the logical OR between two predicates x, y
*      Assumes three computing parties holding the arithmetic shares:
*      x = x1 + x2 + x3
*      y = y1 + y2 + y3
*
*      P1 has: x1, x2, y1, y2
*      P2 has: x2, x3, y2, y3
*      P3 has: x3, x1, y3, y1
*
*      Predicate x is TRUE if x == 0, FALSE otherwise (same for y)
*      x OR y <=> x * y == 0
*
*      The three parties are positioned on a ring. Each party has generated
*      a seed and sent it to its successor party. Each party ends up with
*      two seeds: its own and the seed from its predecessor on the ring.
*      Each party uses these two seeds to generate two random numbers. Then,
*      each party uses the random numbers to generate its rnum as follows:
*      P1: P1_rnum = random(P1_seed) - random(P3_seed)
*      P2: P2_rnum = random(P2_seed) - random(P1_seed)
*      P3: P3_rnum = random(P3_seed) - random(P2_seed)
*      The following holds: P1_rnum + P2_rnum + P3_rnum = 0
*
*      P1 calls or_a() with input: x1, x2, y1, y2, P1_rnum
*      P2 calls or_a() with input: x2, x3, y2, y3, P2_rnum
*      P3 calls or_a() with input: x3, x1, y3, y1, P3_rnum
*
*      Does NOT require Beaver triples or_a all-to-all communication
*      BUT requires pseudo-random sharing of zero
* *****************************************************************/
AShare or_a(AShare x1, AShare x2, AShare y1, AShare y2, AShare rnum);

/*******************************************************************
* and: Computes the logical AND between two predicates x, y
*      Assumes three computing parties holding the arithmetic shares:
*      x = x1 + x2 + x3
*      y = y1 + y2 + y3
*
*      P1 has: x1, x2, y1, y2
*      P2 has: x2, x3, y2, y3
*      P3 has: x3, x1, y3, y1
*
*      Predicate x is TRUE if x == 0, FALSE otherwise (same for y)
*      x AND y <=> x^2 + y^2 == 0
*
*      The three parties are positioned on a ring. Each party has generated
*      a seed and sent it to its successor party. Each party ends up with
*      two seeds: its own and the seed from its predecessor on the ring.
*      Each party uses these two seeds to generate two random numbers. Then,
*      each party uses the random numbers to generate its rnums as follows:
*      P1: P1_rnum1 = random(P1_seed) - random(P3_seed)
*      P1: P1_rnum2 = random(P1_seed) - random(P3_seed) (2nd call to random)
*      P2: P2_rnum1 = random(P2_seed) - random(P1_seed)
*      P2: P2_rnum2 = random(P2_seed) - random(P1_seed) (2nd call to random)
*      P3: P3_rnum = random(P3_seed) - random(P2_seed)
*      P3: P3_rnum2 = random(P3_seed) - random(P2_seed) (2nd call to random)
*      The following hold: P1_rnum1 + P2_rnum1 + P3_rnum1 = 0
*                          P1_rnum2 + P2_rnum2 + P3_rnum2 = 0
*
*      P1 calls and() with input: x1, x2, y1, y2, P1_rnum1, P1_rnum2
*      P2 calls and() with input: x2, x3, y2, y3, P2_rnum, P2_rnum2
*      P3 calls and() with input: x3, x1, y3, y1, P3_rnum, P3_rnum2
*
*      Does NOT require Beaver triples or all-to-all communication
*      BUT requires pseudo-random sharing of zero
*
* *****************************************************************/
AShare and_a(AShare x1, AShare x2, AShare y1, AShare y2, AShare rnum1,
             AShare rnum2);

/*******************************************************************
 * not_b: Computes the complement of a boolean share x
 *      Result: NOT x
 *      Does not_b require party-to-party communication
 * *****************************************************************/
BShare not_b(BShare x);

// TODO: inline (and_b_all_group - and_b - ltz_b - and_b_all)

/*******************************************************************
* and_b: Computes the logical AND between two numbers x, y
*      Assumes three computing parties holding the boolean shares:
*      x = x1 ^ x2 ^ x3
*      y = y1 ^ y2 ^ y3
*
*      P1 has: x1, x2, y1, y2
*      P2 has: x2, x3, y2, y3
*      P3 has: x3, x1, y3, y1
*
*      The three parties are positioned on a ring. Each party has generated
*      a seed and sent it to its successor party. Each party ends up with
*      two seeds: its own and the seed from its predecessor on the ring.
*      Each party uses these two seeds to generate two random numbers. Then,
*      each party uses the random numbers to generate its rnums as follows:
*      P1: P1_rnum = random(P1_seed) ^ random(P3_seed)
*      P2: P2_rnum = random(P2_seed) ^ random(P1_seed)
*      P3: P3_rnum = random(P3_seed) ^ random(P2_seed)
*      The following holds: P1_rnum ^ P2_rnum ^ P3_rnum = 0
*
*      P1 calls and() with input: x1, x2, y1, y2, P1_rnum
*      P2 calls and() with input: x2, x3, y2, y3, P2_rnum
*      P3 calls and() with input: x3, x1, y3, y1, P3_rnum
*
*      Does NOT require Beaver triples or all-to-all communication
*      BUT requires pseudo-random sharing of zero
*
* *****************************************************************/
BShare and_b(const BShare x1, const BShare x2, const BShare y1,
             const BShare y2, const BShare rnum);

/*******************************************************************
* and_b_array: Computes the logical AND between two arrays of
        boolean shares and returns the result in res.

* *****************************************************************/
void and_b_array(const BShare *x1, const BShare *x2,
                 const BShare *y1, const BShare *y2,
                 const BShare *rnum, int len, BShare *res);

/*******************************************************************
* and_b_all: Computes the logical AND between the elements of an array
        in logN communication rounds.
*
*       Modifies input arrays in-place and returns the result shares
*       as the first elements of x1, x2
* *****************************************************************/
void and_b_all(BShare *x1, BShare *x2, int len);

/*******************************************************************
* and_b_all_group: Computes the logical AND between the elements of 
*                  groups in an array in logN communication rounds,
*                  where N is the group size
*
*       - num_groups: the total number of groups in the array
*       - group_size: the number of elements to AND per group
*
*       Modifies input array in-place and returns the result shares
*       as the first 'num_groups' elements of x1, x2
* *****************************************************************/
void and_b_all_group(BShare *x1, BShare *x2, int num_groups, int group_size);

/*******************************************************************
* and_bit: Same as and_b() above but works with single-bit shares
*
* *****************************************************************/
BitShare and_bit(BitShare x1, BitShare x2, BitShare y1, BitShare y2,
                 BitShare rnum);

/*******************************************************************
* and_bit_array: and_bit's corresponding array-based method.

* *****************************************************************/
void and_bit_array(const BitShare *x1, const BitShare *x2,
                        const BitShare *y1, const BitShare *y2,
                        const BitShare *rnum, int len, BitShare *res);

/*******************************************************************
* cmp_swap: Compares and swaps two numbers x, y using their boolean shares
*           Returns:
*             - x = min{x,y}
*             - y = max{x,y}
*
*           Requires computing the following formulas:
*             - min = b * y + (1-b) * x
*             - max = b * x + (1-b) * y
*             where b = x ?> y
*
*           rnums is an array of 319 + 4 = 323 random numbers
*
*           Requires 3 + logN rounds of communication, where N is the length
*           of x, y in number of bits
* *****************************************************************/
void cmp_swap(BShare* x1, BShare* x2, BShare* y1, BShare* y2,
              const BShare* rnums);

/*******************************************************************
* cmp_swap_g: Generalized cmp_swap()
*           Expects two arrays of BShares, compares the elements and swaps the
*           arrays accordingly
*
*           Assumes that the input arrays contain the two shares for each secret
*           at positions att_index1, att_index1+1 and att_index2, att_index2+1
*
*           Returns:
*             - array_1 = the array with the min element
*             - array_2 = the array with the max element
*
*           Requires computing the following formulas:
*             - min = b * r2 + (1-b) * r1
*             - max = b * r1 + (1-b) * r2
*             where b = r1[x] ?> r2[y]
*
*           Requires 3 + logN rounds of communication, where N is the length
*           of each element in r1, r2 in number of bits
* *****************************************************************/
void cmp_swap_g(BShare* r1, BShare* r2, unsigned att_index1,
                unsigned att_index2, int num_elements, int asc);

/*******************************************************************
 * xor_b: Computes the logical XOR between two boolean shares x, y
 *      Result: x ^ y
 *      Does not require party-to-party communication
 * *****************************************************************/
BShare xor_b(BShare x, BShare y);

/*******************************************************************
 * geq: Checks x >= y using boolean shares, where x,y have the same sign
 *      Assumes three computing parties holding the boolean shares:
 *      x = x1 ^ x2 ^ x3
 *      y = y1 ^ y2 ^ y3
 *
 *      P1 has: x1, x2, y1, y2
 *      P2 has: x2, x3, y2, y3
 *      P3 has: x3, x1, y3, y1
 *
 *      x >= y <=> (x_l ^ y_l) & x_l
 *      ^ ~(x_l ^ y_l) & (x_{l−1} ^ y_{l−1}) & x_{l−1}
 *      ^ ~(x_l ^ y_l) & ~(x_{l−1} ^ y_{l−1}) & (x_{l−2} ^ y_{l−2}) & x_{l−2}
 *      ^ ...
 *      ^ ~(x_l ^ y_l) & ~(x_{l−1} ^ y_{l−1}) &...& ~(x_2 ^ y_2) & ~(~x_1 & y_1)
 *
 *      The three parties are positioned on a ring. Each party has generated
 *      a seed and sent it to its successor party. Each party ends up with
 *      two seeds: its own and the seed from its predecessor on the ring.
 *      Each party uses these two seeds to generate random numbers. Then,
 *      each party uses the random numbers to generate each rnum as follows:
 *      P1: P1_rnum = random(P1_seed) ^ random(P3_seed)
 *      P2: P2_rnum = random(P2_seed) ^ random(P1_seed)
 *      P3: P3_rnum = random(P3_seed) ^ random(P2_seed)
 *      The following holds: P1_rnum ^ P2_rnum ^ P3_rnum = 0
 *
 *
 *      P1 calls and() with input: x1, x2, y1, y2, P1_rnums
 *      P2 calls and() with input: x2, x3, y2, y3, P2_rnums
 *      P3 calls and() with input: x3, x1, y3, y1, P3_rnums
 *
 *      Does NOT require Beaver triples or all-to-all communication
 *      BUT requires pseudo-random sharing of zero
 *
 *      Requires 1+log(N) rounds of communication, i.e 7 rounds for 64-bit ints
 *
 * *****************************************************************/
BitShare geq(const BShare x1, const BShare x2,
             const BShare y1, const BShare y2);

/*******************************************************************
* geq_batch: Same as geq() but works in batch mode
*
* *****************************************************************/
void geq_batch(const BShare *x1, const BShare *x2, const BShare *y1,
               const BShare *y2, int numElements, BitShare *res);

/*******************************************************************
 * greater: Checks x > y using boolean shares, where x,y have the same sign
 *      Assumes three computing parties holding the boolean shares:
 *      x = x1 ^ x2 ^ x3
 *      y = y1 ^ y2 ^ y3
 *
 *      P1 has: x1, x2, y1, y2
 *      P2 has: x2, x3, y2, y3
 *      P3 has: x3, x1, y3, y1
 *
 *      x > y <=> (x_l ^ y_l) & x_l
 *       ^ ~(x_l ^ y_l) & (x_{l−1} ^ y_{l−1}) & x_{l−1}
 *       ^ ~(x_l ^ y_l) & ~(x_{l−1} ^ y_{l−1}) & (x_{l−2} ^ y_{l−2}) & x_{l−2}
 *       ^ ...
 *       ^ ~(x_l ^ y_l) & ~(x_{l−1} ^ y_{l−1}) &...& ~(x_2 ^ y_2) & (x_1 & ~y_1)
 *
 *      The three parties are positioned on a ring. Each party has generated
 *      a seed and sent it to its successor party. Each party ends up with
 *      two seeds: its own and the seed from its predecessor on the ring.
 *      Each party uses these two seeds to generate random numbers. Then,
 *      each party uses the random numbers to generate each rnum as follows:
 *      P1: P1_rnum = random(P1_seed) ^ random(P3_seed)
 *      P2: P2_rnum = random(P2_seed) ^ random(P1_seed)
 *      P3: P3_rnum = random(P3_seed) ^ random(P2_seed)
 *      The following holds: P1_rnum ^ P2_rnum ^ P3_rnum = 0
 *
 *      The rnums array stores all rnums we need for evaluating the above
 *      formula. For 64-bit integers, we need 128+191=319 rnums per party.
 *
 *      P1 calls and() with input: x1, x2, y1, y2, P1_rnums
 *      P2 calls and() with input: x2, x3, y2, y3, P2_rnums
 *      P3 calls and() with input: x3, x1, y3, y1, P3_rnums
 *
 *      Does NOT require Beaver triples or all-to-all communication
 *      BUT requires pseudo-random sharing of zero
 *
 *      Requires 1+log(N) rounds of communication, i.e 7 rounds for 64-bit ints
 * *****************************************************************/
BitShare greater(const BShare x1, const BShare x2,
                 const BShare y1, const BShare y2);

/*******************************************************************
* greater_batch: Same as greater() but works in batch mode
*
* *****************************************************************/
void greater_batch(const BShare *x1, const BShare *x2,
                   const BShare *y1, const BShare *y2,
                   int numElements, BitShare *res);

/*******************************************************************
* greater_batch2: Same as greater_batch() but compares elements in a
*       single array
*
*       Used by group_by_count_sel_odd_even()
* *****************************************************************/
void greater_batch2(BShare **c, int idx, int start, int numElements, 
                    int dist, BitShare *res);

/****************************************************************************
 * eq_b: Checks x == y using boolean shares.
 *      Assumes three computing parties holding the boolean shares:
 *      x = x1 ^ x2 ^ x3
 *      y = y1 ^ y2 ^ y3
 *
 *      P1 has: x1, x2, y1, y2
 *      P2 has: x2, x3, y2, y3
 *      P3 has: x3, x1, y3, y1
 *
 *      The equality returns 1 if all corresponding bits in the inputs
 *      are the same and 0 otherwise, using this formula:
 *      x == y <=> (x_l ^ y_l ^ 1) & (x_{l-1} ^ y_{l-1} ^ 1) & ...
 *                  & (x_1 ^ y_1 ^ 1)
 *
 *      The computation requires 6 rounds of communication.
 *      Each party initially computes res0_i = x_i ^ y_i ^ 1 locally
 *      for i=0 to 63 (for 64-bit ints) and for both shares of each input.
 *      It then executes the and operations in a tree with log2(64) = 6 levels:
 *      - At L=1: 32 concurrent and operations where res1_i = res0_i ^ res0_i+1
 *      - At L=2: 16 concurrent and operations where res2_i = res1_i ^ res1_i+1
 *      - At L=3: 8 concurrent and operations where res3_i = res2_i ^ res2_i+1
 *      - At L=4: 4 concurrent and operations where res4_i = res3_i ^ res3_i+1
 *      - At L=5: 2 concurrent and operations where res5_i = res4_i ^ res4_i+1
 *      - At L=6: 1 and operation where res6_i = res5_i ^ res5_i+1
 *
 *      The three parties are positioned on a ring. Each party has generated
 *      a seed and sent it to its successor party. Each party ends up with
 *      two seeds: its own and the seed from its predecessor on the ring.
 *      Each party uses these two seeds to generate random numbers. Then,
 *      each party uses the random numbers to generate each rnum as follows:
 *      P1: P1_rnum = random(P1_seed) ^ random(P3_seed)
 *      P2: P2_rnum = random(P2_seed) ^ random(P1_seed)
 *      P3: P3_rnum = random(P3_seed) ^ random(P2_seed)
 *      The following holds: P1_rnum ^ P2_rnum ^ P3_rnum = 0
 *
 *      The rnums array stores all rnums we need for evaluating the above
 *      formula. For 64-bit integers, we need 64-1 = 63 rnums per party.
 *
 *      P1 calls and() with input: x1, x2, y1, y2, P1_rnums
 *      P2 calls and() with input: x2, x3, y2, y3, P2_rnums
 *      P3 calls and() with input: x3, x1, y3, y1, P3_rnums
 *
 *      Does NOT require Beaver triples or all-to-all communication
 *      BUT requires pseudo-random sharing of zero
 *
 *      Requires logN rounds of communication, i.e., 6 rounds for 64-bit ints
 * **************************************************************************/
BShare eq_b(BShare x1, BShare x2, BShare y1, BShare y2);

// synchronous element-wise version
BShare eq_b_sync(BShare x1, BShare x2, BShare y1, BShare y2);

// asynchronous element-wise
BShare eq_b_async(BShare x1, BShare x2, BShare y1, BShare y2);

/****************************************************************************
 * eq_b_array: Array-based implementation of eq_b().
 * **************************************************************************/
void eq_b_array(BShare *x1, BShare *x2, BShare *y1, BShare *y2, long len, BShare *res);

// array-based equality with interleave between levels
void eq_b_array_inter(BShare *x1, BShare *x2, BShare *y1, BShare *y2, long len, BShare *res);

// array-based boolean equality with asynchronous exchange
// and interleaving between batches of elements
void eq_b_array_inter_batch(BShare *x1, BShare *x2, BShare *y1, BShare *y2, long len,
                BShare *res);

/****************************************************************************
* Computes bitwise equality for one tree level.
* ****************************************************************************/
BShare eq_b_level2(int numbits, BShare z1, BShare z2);

/*******************************************************************
* boolean_addition: Adds two numbers x and y using their boolean shares and
*                   and stores the boolean shares of z = x + y in z1 and z2
*
*      The boolean addition of two numbers is computed as follows:
*
*      1. Perform a bitwise XOR:
*         - z1 = x1 ^ y1
*         - z2 = x2 ^ y2
*
*      2. Compute the carry at position i, 0 <= i < 64 using the recursive
*         formula:
*         - c_0 = 0  (carry at position 0)
*         - c_1 = x1 AND y1  (carry at position 1)
*         - c_i = (x_i AND y_i) ^ (x_i AND c_{i-1}) ^ (y_i AND c_{i-1}), i > 1
*
*      3. Do a final XOR:
*         - z1 ^= c
*         - z2 ^= c
*
*      This function assumes three computing parties holding the shares:
*      x = x1 ^ x2 ^ x3
*      y = y1 ^ y2 ^ y3
*
*      P1 has: x1, x2, y1, y2
*      P2 has: x2, x3, y2, y3
*      P3 has: x3, x1, y3, y1
*
*      rnums is an array of 187 (= 1 + 3*62) random numbers needed for the
*      boolean addition, which are generated by each party as follows:
*
*      The three parties are positioned on a ring. Each party has generated
*      a seed and sent it to its successor party. Each party ends up with
*      two seeds: its own and the seed from its predecessor on the ring.
*      Each party uses these two seeds to generate random numbers. Then,
*      each party uses the random numbers to generate each rnum as follows:
*      P1: P1_rnum = random(P1_seed) ^ random(P3_seed)
*      P2: P2_rnum = random(P2_seed) ^ random(P1_seed)
*      P3: P3_rnum = random(P3_seed) ^ random(P2_seed)
*      The following holds: P1_rnum ^ P2_rnum ^ P3_rnum = 0
*
*
*      The function requires N rounds of communication, where N is the
*      length of each boolean share in number of bits. At each round, we need
*      to exchange two BShares between each pair of parties
* *****************************************************************/
void boolean_addition(BShare x1, BShare x2, BShare y1, BShare y2,
                      BShare *z1, BShare *z2,
                      BShare *rnums);

/****************************************************************************
 * ltz_b: Checks x<0 using boolean shares.
 *      Assumes three computing parties holding the boolean shares:
 *      x = x1 ^ x2 ^ x3
 *      y = y1 ^ y2 ^ y3
 *
 *      P1 has: x1, x2, y1, y2
 *      P2 has: x2, x3, y2, y3
 *      P3 has: x3, x1, y3, y1
 *
 *      This is a local operation that simply checks whether the MSB is 1.
 *
 * **************************************************************************/
BShare ltz_b(BShare x);

/****************************************************************************
 * ltz_b_array: Checks x<0 using boolean shares for an array of elements.
 *
 * **************************************************************************/
void ltz_b_array(const BShare *x, int len, BShare *res);

/****************************************************************************
* boolean_addition_batch: Batch version of boolean_addition().
*               Returns a single share of the result.
* **************************************************************************/
void boolean_addition_batch(BShare *x1, BShare *x2, BShare *y1, BShare *y2,
                       BShare *res1, int numElements);

/****************************************************************************
* boolean_addition_batch2: Batch version of boolean_addition().
*       Returns a single share of the result in res1.
*
*       Used by group_by_sum_rca_sel_odd_even() and group_by_sum_rca_odd_even()
* **************************************************************************/
void boolean_addition_batch2(BShare **c, BShare *res1, 
                             int start, int dist, int numElements,
                             int idx);

/****************************************************************************
 * convert_single_bit: Converts a single-bit BShare to an AShare.
 *
 * NOTE: During pre-processing, the parties must have received a set of random
 *       bit shares generated by generate_rand_bit_shares().
 *
 *      The method follows the corresponding Crypten implementation
 *      https://github.com/facebookresearch/CrypTen
 *      and works as follows:
 *
 *      1. Each party computes [z] = bit^[rb] and reveals the value z.
 *         To do that, P2 and P3 send their shares [z] to P1, and P1 sends back
 *         the result z.
 *      2. The arithmetic share is [ra] = [ra] * (1-2*z) + z.
 * **************************************************************************/
AShare convert_single_bit(BShare bit, AShare ra, BShare rb);


/****************************************************************************
 * convert_single_bit_array: Array-based implementation of convert_single_bit()
 *         The result is stored in res.
 * **************************************************************************/
void convert_single_bit_array(BShare *bit, AShare *ra, BShare *rb, int len,
                               AShare *res);

/****************************************************************************
 * convert_a_to_b_array: Converts an array of arithmetic shares to boolean
 *                      shares using the protocol described in ABY^3
 *                      https://eprint.iacr.org/2018/403.pdf.
 *
 *      Consider [[x]]a = (x1, x2, x3) where x = x1+x2+x3.
 *      Since we use replicated sharing, P1 holds both x1 and x2
 *      and can compute x1 + x2 locally.
 *
 *      P1 generates boolean shares of (x1+x2) and P3 generates boolean
 *      shares of x3 and send parties their shares.
 *
 *      We then use the ripple-carry-adder to compute:
 *      [[rca_res]]_b = ([[x1+x2]]_b + [[x3]]_b).
 *
 *      Each party also generates two random shares [[r1]]_b and [[r2]]_b
 *      and reveals them to parties (1, 2) and (2, 3) respectively.
 *
 *      Then, each party locally computes
 *      [[y1]]_b = [[rca_res]]_b ^ [[r1]]_b ^ [[r2]]_b and reveals the result
 *      to parties 1, 3.
 *
 *      [[x]]_b = (y1, y2, y3)
 *
 *      - xa1: array of local arithmetic shares
 *      - xa2: array of remote arithmetic shares
 *      - xb1: local share of result
 *      - xb2: remote share of result
 *      - len: array length
 * **************************************************************************/
void convert_a_to_b_array(AShare *xa1, AShare *xa2, BShare *xb1, BShare *xb2,
                          int len);

/****************************************************************************
* cmp_swap_batch: Array-based implementation of cmp_swap_g()
*         Updates intput table in place
*
*         Returns the index of the last row checked in the batch
*
*         - first_index: start of the current batch in 'rows'
*         - rows: the rows (content) of the shares table
*         - length: the total number of rows in the shares table
*         - row_length: the row length (incl. remote shares of each attribute)
*         - att_indices: the indices of the attribute to sort by
*         - num_attributes: the size of 'att_indices' array
*         - phase: the phase of the bitonic sort. There are logN phases in
*                  total, where N is the number of rows (i.e. 'length')
*         - column: the column within each phase. Phase 0 has 1 column, Phase 2
*                   has 2 columnes, ..., phase logN has longN columns
*         - asc: the sorting direction (ASC/DESC) of each attribute (ASC=1)
*         - num_comparisons: the batch size (<= length / 2)
*
*         The number of rounds depends on the number of attributes:
*         - num_attributes=1: log(L+1) + 1
*         - num_attributes=2: num_attributes * (log(L+1) * logL) + 1 + 2
*         - num_attributes=3: num_attributes * (log(L+1) * logL) + 1 + 4
*
*         where L is the length of the sorting attributes in number of bits
* **************************************************************************/
int cmp_swap_batch(int first_index, BShare** rows, int length, int row_length,
                   unsigned* att_indices, int num_attributes, int phase,
                   int column, bool* asc, int num_comparisons);

#endif
