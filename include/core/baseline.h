#ifndef BASELINE_H
#define BASELINE_H

#include <stdio.h>
#include <string.h>
#include <math.h>

#include "comm.h"
#include "mpctypes.h"
#include "party.h"
#include "utils.h"


/*******************************************************************
 * group_by_sum_rca: group-by-sum that uses a Ripple-Carry Adder
 *          Expects a BShareTable sorted on the group-by key(s).
 *
 *          The two boolean shares of the sum attribute are expected to be in
 *          the last two columns of the given table.
 *
 *          The function updates the sum attribute in place and does not account
 *          for 'selected' bits.
 *
 *          - key_indices: An array of the group-by keys' indices
 *          - num_keys: The number of group-by keys
 * *****************************************************************/
void group_by_sum_rca(BShareTable* table, unsigned* key_indices,
                                          unsigned num_keys);

/*******************************************************************
 * group_by_sum_rca_sel: Same as group_by_sum_rca() but takes into account the
 *          selected bits. Expects a BShareTable sorted on the group-by key(s).
 *
 *          The two boolean shares of the sum attribute are expected to be in
 *          the last two columns of the given table.
 *
 *          - selected: The array of 'selected' bits
 *          - key_indices: An array of the group-by keys' indices
 *          - num_keys: The number of group-by keys
 * *****************************************************************/
void group_by_sum_rca_sel(BShareTable* table, BShare* selected,
                          unsigned* key_indices, unsigned num_keys);

/*******************************************************************
 * group_by_sum_rca_sel_odd_even: Same as group_by_sum_rca_sel() but uses
 *          odd-even aggregation to reduce the number of rounds.
 *
 *          The two boolean shares of the sum attribute are expected to be in
 *          the last two columns of the given table.
 *
 *          - batch_size: The number of elements exchanged per round
 *          - selected_b: The array of 'selected' bits
 *          - key_indices: An array of the group-by keys' indices
 *          - num_keys: The number of group-by keys
 * *****************************************************************/
void group_by_sum_rca_sel_odd_even(BShareTable* table, int batch_size, 
                                   BShare* selected_b,
                                   unsigned* key_indices, unsigned num_keys);

/*******************************************************************
 * group_by_sum_rca_sel_odd_even: Same as group_by_sum_rca_sel() but uses
 *          odd-even aggregation to reduce the number of rounds.
 *
 *          The two boolean shares of the sum attribute are expected to be in
 *          the last two columns of the given table.
 *
 *          The function updates the sum attribute in place and does not account
 *          for 'selected' bits.
 * 
 *          - key_indices: An array of the group-by keys' indices
 *          - num_keys: The number of group-by keys
 * *****************************************************************/
void group_by_sum_rca_odd_even(BShareTable* table,  int batch_size, 
                               unsigned* key_indices, unsigned num_keys);

/*******************************************************************
* select_greater_batch_const: Compares all elements in the specified column
*          with the given secret-shared constant.
*
*          Used in password reuse and credit score queries
*
*          Returns the result bits in result.
* *****************************************************************/
void select_greater_batch_const(BShareTable input, int leftcol,
                                BShare cshare1, BShare cshare2,
                                BShare result[]);

/*******************************************************************
* select_eq_batch_const: Compares all elements in the specified column
*          with the given secret-shared constant.
*
*          Used in credit score query (baseline)
*
*          Returns the result bits in result.
* *****************************************************************/
void select_eq_batch_const(BShareTable input, int leftcol,
                           BShare cshare1, BShare cshare2,
                           BShare result[]);

/****************************************************************************
* boolean_addition_batch_const: Adds all elements in the specified column with
*               the given secret-shared constant.
* **************************************************************************/
void boolean_addition_batch_const(BShareTable input, int leftcol,
                      BShare cshare1, BShare cshare2, BShare *res);

/*******************************************************************
 * geq_batch_const: geq_batch() with a public constant
 *          Compares all elements in the given array with the given
 *          secret-shared constant.
 *
 *          Used in test_tpch_q4_baseline.
 *
 *          Returns the result bits in res.
 * *****************************************************************/
void geq_batch_const(const BShare *x1, const BShare *x2, BShare y1, BShare y2,
                     int numElements, BitShare *res);

/*******************************************************************
* greater_batch_const: greater_batch() with a public constant
*          Compares all elements in the given array with the given
*          secret-shared constant.
*
*          Used in test_tpch_q4_baseline.
*
*          Returns the result bits in res.
* *****************************************************************/
void greater_batch_const(const BShare *x1, const BShare *x2, BShare y1,
                         BShare y2, int numElements, BitShare *res);

/*******************************************************************
 * eq_b_array_cont: boolean equality between a vector and a
 *                  secret-shared constant.
 *
 *          Used in test_tpch_q13_baseline.
 *
 *          Returns the result bits in res.
 * *****************************************************************/
void eq_b_array_const(BShare *x1, BShare *x2, BShare y1, BShare y2, long len,
                BShare *res);

#endif
