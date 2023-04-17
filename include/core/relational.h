#ifndef RELATIONAL_H
#define RELATIONAL_H

#include "assert.h"
#include "mpctypes.h"
#include "utils.h"
#include "party.h"
#include "primitives.h"
#include "baseline.h"
#include <stdio.h>

/**
 * An enumeration of the supported operations inside predicates.
 * EQ: equality
 * LT: less than
 * GT: greater than
 * LEQ: less than or equal
 * GEQ: greater than or equal
 * EQC: equality with public constant
 * GEQC: grearter than or equal to public constant
 * GC: greater than public constant
 * GR: greater than between two columns
 * EXP: expression
 **/
typedef enum {EQ, LT, GT, LEQ, GEQ, EQC, GEQC, GC, GR, EXP} OP_TYPE;

/**
 * Predicate represents an equality, inequality, or expression
 * predicate provided to selection and join operations.
 * A predicate is composed by two operands, op1 and op2,
 * which can themselves be predicates, and an operation type.
 * When op1 and/or op2 are NULL, then left and right point to
 * the input shares to be used for this op_type.
**/
struct predicate {
    OP_TYPE operation;
    struct predicate *op1;
    struct predicate *op2;
    int leftcol;
    int rightcol;
};

/**
 * The corresponding predicate for boolean sharing.
**/
struct predicate_b {
    OP_TYPE operation;
    struct predicate_b *op1;
    struct predicate_b *op2;
    int leftcol;
    int rightcol;
    BShare cs1;   // Used only by EQC, GEQC, and GC
    BShare cs2;   // Used only by EQC, GEQC, and GC
};

typedef struct predicate Predicate;
typedef struct predicate_b Predicate_B;

/*******************************************************************
 * select_a: Performs a selection to the input table share using the
 *         predicate p. Each input row produces an output which is
 *         stored in the result array.
 * *****************************************************************/
void select_a(AShareTable input, Predicate p, AShare c1, AShare c2, AShare result[]);

void and_b_table(const BShareTable input, int leftcol, int rightcol, int size, BShare result[]);

/*******************************************************************
 * select_b: Performs a selection on the input table share using the
 *           predicate p. The result array contains one value (0/1)
 *           per input row, indicating whether the corresponding row
 *           should be included in the result.
 *           This method does not perform a projection.
 *
 * We currently support equality predicates, e.g.:
 *
 *      SELECT id FROM r1 WHERE r1.att = c
 *
 * where c is a public constant.
 *
 * To compute equality with a constant, we use the following protocol:
 *   1. During pre-processing, the data owner sends shares of
 *      `att-c` and `c-att` to the computing parties.
 *   2. The parties compute `z1 = (c-att < 0)` and `z2 = (att-c < 0)`
 *      using a local call of the ltz_b() primitive.
 *   3. The result is r = z1 ^ z2 ^ 1:
 *      - z1 ^ z2 = 0, if att=c
 *      - z1 ^ z2 = 1, if either c-att < 0 or att-c < 0.
 *
 * `c-att` and `att-c` need to be included in the input BShareTable
 *  at the columns pointed to by the predicate's leftcol and rightcol
 *  fields.
 *
 *
 * *****************************************************************/
void select_b(BShareTable input, Predicate_B p, BShare result[]);

/*******************************************************************
 * join_b: Performs a join on the input table shares using the
 *         predicate p. Each input row produces an output which is
 *         stored in the result array.
 * *****************************************************************/
void join_b(BShareTable input1, BShareTable input2, Predicate_B p, BShare result[]);

/*******************************************************************
 * join_b_batch: Performs a batched join on the input table shares
 *          using the provided table indices indicating the batch
 *          start and end locations.
 *          The predicate contains the column indexes of the join
 *          attributes.
 * *****************************************************************/
void join_b_batch(BShareTable *input1, BShareTable *input2,
                    int start1, int end1, int start2, int end2,
                    Predicate_B p, BShare *remote, BShare *result);

void join_eq_b_batch(BShareTable *input1, BShareTable *input2,
                     int start1, int end1, int start2, int end2,
                     int leftcol, int rightcol, BShare* remote,
                     BShare* result);

void join_geq_b_batch(BShareTable *input1, BShareTable *input2,
                              int start1, int end1, int start2, int end2,
                              int leftcol, int rightcol, BShare* _remote,
                              BShare* result);

/*******************************************************************
* distinct: Computes distinct on the given BShareTable.
*
*           Expects a sorted BShareTable.
*
*           - att_index is the index of the DISTINCT attribute
*           - distinct is an array of bit shares indicating the rows with
*             distinct elements
*
*           Relies on eq_b() and requires longN communication rounds,
*           where N is the length of the elements compared with eq_b()
*           in number of bits
*
* *****************************************************************/
void distinct(BShareTable* table, unsigned att_index, BitShare* distinct);


/*******************************************************************
* distinct_batch: Same as distinct but works in batch mode
*
*           Expects a sorted BShareTable.
*
*           - key_indices are the indixes of the DISTINCT keys
*           - num_keys is the number of DISTINCT keys
*           - distinct is an array of bit shares indicating the rows with
*             distinct elements
*           - num_comparisons is the number of comparisons per batch
*
*           Relies on eq_b() and requires longN communication rounds,
*           where N is the length of the elements compared with eq_b()
*           in number of bits
*
* *****************************************************************/
void distinct_batch(BShareTable* table, unsigned* key_indices, unsigned num_keys,
                    BitShare* distinct, unsigned num_comparisons);

/*******************************************************************
* distinct_linear: Computes distinct on the given BShareTable
*                   taking previous selection into account
*                   with a linear scan.
*
*           Expects a sorted BShareTable on the DISTINCT key(s).
*           Ovverides the selected_b column with the resulting distinct bit.
*
*           - key_indices includes the indices of the distinct keys
*           - num_keys is the length of key_indices (number of distinct keys)
*           - selected_b indicates which rows have been selected from previous
*             steps (BShares). The result (distinct bit) is stored in this array.
*
*           Requires O(N) communication rounds, where N is
*           the number of rows in the given table
* *****************************************************************/
void distinct_linear(BShareTable* table, unsigned* key_indices, int num_keys,
                    BShare* selected_b);

// TODO: make this private
void eq_bulk(int, int, BShare*, BShare*, int);

/*******************************************************************
* in: Oblivious semi-join between left and right
*
*           Populates array 'in' with bitshares, indicating the attributes at
*           'att_index' of the left row appears in the right input
*
*           - left is the left relation
*           - right is the right relations (the smaller one)
*             Currently, the size of the right input MUST be a power of two
*           - left_index is the index of the join attribute in the left input
*           - right_index is the index of the join attribute in the right input
*           - in is a bit shares array with size equal to the size of left input
*           - num_rows_left is the number of rows from the left input per batch
*
*           Total number of comparisons per batch is:
*             - First round: O(num_rows_left) * right->numRows +
*             - Second round: O(num_rows_left) * (right->numRows - 1)
*
*           Uses on join_eq_b_batch() and requires num_rows_left * (logL + LogR)
*           communication rounds, where L is the length of the join elements
*           in number of bits, and R is the number of rows in the right input
*
* *****************************************************************/
void in(BShareTable* left, BShareTable* right, unsigned left_index,
        unsigned right_index, BShare* in, unsigned num_rows_left);

/*******************************************************************
* in_sel_right: Fused semi-join selection on the right input
*           - sel_index is the selected bit index in right
*
* *****************************************************************/
void in_sel_right(BShareTable* left, BShareTable* right,
        unsigned left_index, unsigned right_index, unsigned sel_index,
        BShare* res_in, unsigned num_rows_left);

// TODO: group by should affect the [SEL] not the grouping attributes values
//  - Otherwise, these values should be non negative to avoid -1 value confussion
//  with no selection.

/*******************************************************************
* group_by_count: Groups rows of the given table and counts the rows per group
*
*           Expects a sorted BShareTable on the group-by key(s).
*           Updates table in place.
*
*           - key_indices includes the indices of the group-by keys
*           - num_keys is the length of key_indices (number of group-by keys)
*           - selected_b indicates which rows have been selected from previous
*             steps (BShares)
*           - selected indicates which rows have been selected from previous
*             steps (AShares)
*           - rb and ra are arrays of random binary and arithemtic bit shares of
*             length 2*(N-1), where N is the number of rows in the given table
*             these random numbers are needed to compute arithmetic shares from
*             boolean shares
*
*           Requires 2 + N*(6 + logL) = O(N) communication rounds, where N is
*           the number of rows in the given table and L is the length of each
*           table element in number of bits
* *****************************************************************/
void group_by_count(BShareTable* table, unsigned* key_indices, int num_keys,
                    BShare* selected_b, AShare* selected, BShare* rb,
                    AShare* ra);

/*******************************************************************
* group_by_min_max_sel: Groups rows of the given table and computes min and max
*           for the specified attributes. It also takes into account 'selected'
*           bits. Expects a sorted BShareTable on the group-by key(s) and
*           updates the input table in place.
*
*           Used in credit score query
*
*           - selected indicates which rows have been selected from previous
*             steps (BShares)
*           - min_att is the index of the attribute where MIN() is applied
*           - max_att is the index of the attribute where MAX() is applied
*           - key_indices includes the indices of the group-by keys
*           - num_keys is the length of key_indices (number of group-by keys)
*
* *****************************************************************/
void group_by_min_max_sel(BShareTable* table, BShare* selected,
                          unsigned min_att, unsigned max_att,
                          unsigned* key_indices, int num_keys);

/*******************************************************************
* group_by_min_max_sel_odd_even: Groups rows of the given table and computes 
*           min and max for the specified attributes. It also takes into 
*           account 'selected' bits. Expects a sorted BShareTable on the 
*           group-by key(s) and updates the input table in place
*
*           Uses odd-even aggregation to reduce the number of rounds to a
*           logarithmic factor
*
*           Used in credit score query
*
*           - selected_b indicates which rows have been selected from previous
*             steps (BShares)
*           - batch_size is the max number of elements per communication round
*           - min_att is the index of the attribute where MIN() is applied
*           - max_att is the index of the attribute where MAX() is applied
*           - key_indices includes the indices of the group-by keys
*           - num_keys is the length of key_indices (number of group-by keys)
*
* *****************************************************************/
void group_by_min_max_sel_odd_even(BShareTable* table, unsigned batch_size, 
                                   BShare* selected_b,
                                   unsigned min_att, unsigned max_att,
                                   unsigned* key_indices, int num_keys);

/*******************************************************************
 * group_by_count_micro: group-by version for microbenchmarks.
 *          This method doesn't expect selected bits as input.
 *
 *          - key_indices: An array of the group-by keys' indices
 *          - num_keys: The number of group-by keys
 * *****************************************************************/
void group_by_count_micro(BShareTable* table, unsigned* key_indices,
                          int num_keys, AShare* counters,
                          AShare* remote_counters, BShare* rb, AShare* ra);

/*******************************************************************
 * group_by_sum_rca: group-by-sum that uses RCA.
 *          This method doesn't expect selected bits as input,
 *          i.e. it computes the group sizes.
 *          Expects an input BShareTable sorted on the group-by key(s)
 *
 *          - key_indices: An array of the group-by keys' indices
 *          - num_keys: The number of group-by keys
 * *****************************************************************/
void group_by_sum_rca(BShareTable* table, unsigned* key_indices,
                                          unsigned num_keys);

/*******************************************************************
* group_by_join: Applies a fused group-by-join-aggregation (sum) operation
            on two input tables, left and right.
*
*           Expects a sorted left BShareTable on group_att_index.
*
*           - batch_left:  the number of tuples from the left input that form
*                          a batch.
*           - group_att_index:  the index of the attribute to group by and it
*                               refers to the left input.
*           - left_join_index:  the join attribute on the left input
*           - right_join_index: the join attribute on the right input
*           - right_att_index:  the aggregation attribute on the right input
                                NOTE: this must be an arithmetic share
*           - sum_res_index:  The index of the aggregated values in the
                                left input table.
*           - count_res_index: The index of the counts
*       NOTE: does not exchange remote results internally
* *****************************************************************/
void group_by_join(BShareTable* left, BShareTable* right, int start_left,
                   int end_left, int group_att_index, int left_join_index,
                   int right_join_index, int right_att_index,
                   BShare* rb_left, AShare* ra_left, BShare* rb_right,
                   AShare* ra_right, unsigned sum_res_index,
                   unsigned count_res_index);

/*******************************************************************
* group_by_join_first: Applies first phase of the group-by-join-aggregation
*           (sum) operation on two input tables, left and right. 
*
*           - batch_left:  the number of tuples from the left input that form
*                          a batch.
*           - left_join_index:  the join attribute on the left input
*           - right_join_index: the join attribute on the right input
*           - right_att_index:  the aggregation attribute on the right input
                                NOTE: this must be an arithmetic share
*           - sum_res_index:  The index of the aggregated values in the
                                left input table.
* *****************************************************************/
void group_by_join_first(BShareTable* left, BShareTable* right, int start_left,
                   int end_left, int left_join_index,
                   int right_join_index, int right_att_index, BShare* rb_right,
                   AShare* ra_right, unsigned sum_res_index);

/*******************************************************************
* group_by_sum_odd_even: Applies second phase of the group-by-join-aggregation
*           (sum) operation on two input tables, left and right. 
*
*           - batch_size:  the number of tuples from the left input that form
*                          a batch.
*           - right_att_index:  the aggregation attribute on the right input
                                NOTE: this must be an arithmetic share
*           - sum_res_index:  The index of the aggregated values in the
                                left input table.
*           - count_res_index: The index of the counts (not used)
* *****************************************************************/
void group_by_sum_odd_even(BShareTable* left, int batch_size,
                   BShare* rb_left, AShare* ra_left, unsigned sum_res_index,
                   unsigned count_res_index, unsigned* key_indices, unsigned num_keys);

/*******************************************************************
* group_by_avg_odd_even: Computes the average (sum and count) on a sorted relation. 
*
*           Expects the relation to be sorted on group-by-keys (key_indices)
*           Expects column at count_res_index to be initiailized to `1`
*
*           - rel: the input relation
*           - batch_size:  the number of tuples from the left input that form
*                          a batch.
*           - rb_rel: binary shares used in b2a conversion 
*           - ra_rel: arithemic shares used in b2a conversion 
*           - sum_res_index:  The index of the values to sum (in place).
*           - count_res_index: The index to store the counts 
* *****************************************************************/
void group_by_avg_odd_even(BShareTable* rel, int batch_size,
                   BShare* rb_rel, AShare* ra_rel, unsigned sum_res_index,
                   unsigned count_res_index, unsigned* key_indices, unsigned num_keys);

/*******************************************************************
* group_by_count_sel_odd_even: Groups rows of the given table and counts the rows per group
*
*           Uses odd-even aggregation to reduce the number of rounds to a 
*           logarithmic factor. Updates table in place.
*
*           - key_indices includes the indices of the group-by keys
*           - num_keys is the length of key_indices (number of group-by keys)
*           - selected_b indicates which rows have been selected from previous
*             steps (BShares)
*           - selected indicates which rows have been selected from previous
*             steps (AShares)
*           - rb and ra are arrays of random binary and arithemtic bit shares of
*             length 2*(N-1), where N is the number of rows in the given table
*             these random numbers are needed to compute arithmetic shares from
*             boolean shares
*
*           Requires O(nlogn) operations and O(logN) communication rounds, where 
*           N is the number of rows in the given table
* *****************************************************************/
void group_by_count_sel_odd_even(BShareTable* table, unsigned* key_indices, int num_keys,
                        unsigned batch_size, BShare* selected_b, AShare* selected, 
                        BShare* rb, AShare* ra);

/*******************************************************************
* group_by_count_odd_even: Groups rows of the given table and counts the rows per group
*           This method does not expect selected bits.
*
*           Uses odd-even aggregation to reduce the number of rounds to a 
*           logarithmic factor. Updates table in place.
*
*           - key_indices includes the indices of the group-by keys
*           - num_keys is the length of key_indices (number of group-by keys)
*           - rb and ra are arrays of random binary and arithemtic bit shares of
*             length 2*(N-1), where N is the number of rows in the given table
*             these random numbers are needed to compute arithmetic shares from
*             boolean shares
*
*           Requires O(nlogn) operations and O(logN) communication rounds, where 
*           N is the number of rows in the given table
* *****************************************************************/
void group_by_count_odd_even(BShareTable* table, unsigned* key_indices, int num_keys,
                             unsigned batch_size, AShare* counters,
                             AShare* remote_counters,BShare* rb, AShare* ra);

/*******************************************************************
* sort: Sorts a BSharesTable in place using bitonic sort
*
*           Requires a BSharesTable with 2^k number of rows.
*
*           - sort_attribute is the index of the attribute to sort on
*           - desc defines the sorting direction
*           - rnums is an array of 323*N*(logN)^2 random numbers
*
*           Requires (3 + logL)*logN = O(logN) rounds of communication, where L
*           is the length of each table element in number of bits and N is the
*           number of elements (BShares)
* *****************************************************************/
void bitonic_sort(BShareTable* table, int low, int cnt, unsigned sort_attribute,
                  int asc);

/*******************************************************************
* sort: Same as bitonic_sort() but works in batches
*
*       - 'batch_size' and 'table->numRows' must both be a power of 2
*
*       - sort_attributes: the indices of the sort attributes in the BShareTable
*         For example, if we want to sort by the first two attributes, we need
*         to pass an array sort_attributes = {0,2}
*       - num_attributes: the size of the sort_attributes array
*       - asc: The sorting direction (ASC/DESC) for each attribute (ASC=1)
*       - 'batch_size' (<= 'table->numRows'/2) is the number of comparisons
*         within a single round of communication
*
*       Requires logN * (1+ logN)/2 * ceil(N/2*'batch_size') *
                 * rounds_of_cmp_swap_batch (see primitives.h)
*       = O(N * (LogN)^2 / batch_size) communication rounds, where N is the
*       number of rows in the shares table
*
* *****************************************************************/
void bitonic_sort_batch(BShareTable* table, unsigned* sort_attributes,
                        int num_attributes, bool* asc, int batch_size);

/*******************************************************************
* mask: Masks the non-selected rows in a table by updating the table in place.
*
*       - selected contains the selected bit
*       - batch_size defines the number of rows to mask within a single round
* *****************************************************************/
void mask(BShareTable* table, BShare* selected, int batch_size);

void adjacent_geq(BShareTable* table, unsigned att_index1, unsigned att_index2,
                  BitShare* result, unsigned num_comparisons, int swap);

#endif
