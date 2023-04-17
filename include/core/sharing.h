#ifndef SHARING_H
#define SHARING_H

#include <stdbool.h>
#include <sodium.h>
#include "mpctypes.h"
#include "primitives.h"


/*******************************************************************
 * init_sharing: Initializes random seed
 * *****************************************************************/
void init_sharing();

/*******************************************************************
 * generate_int_share: Generates three boolean shares for integer x
 *                            such that x = x1 + x2 + x3
 * *****************************************************************/
void generate_int_share(Data x, AShare *x1, AShare *x2, AShare *x3);

/*******************************************************************
 * generate_bool_share: Generates three boolean shares for integer x
 *                     such that x = x1 ^ x2 ^ x3
 * *****************************************************************/
void generate_bool_share(Data x, BShare *x1, BShare *x2, BShare *x3);

/*******************************************************************
 * generate_bool_share_table: Generates three tables of boolean shares
 *
 *                            Assumes preallocated BShareTables
 * *****************************************************************/
void generate_bool_share_tables(Table *data, BShareTable *shares_1,
                                BShareTable *shares_2, BShareTable *shares_3);

/*******************************************************************
* generate_bool_share_table: Generates three tables of arithmetic shares
*
*                            Assumes preallocated AShareTables
* *****************************************************************/
void generate_int_share_tables(Table *data, AShareTable *shares_1,
                               AShareTable *shares_2, AShareTable *shares_3);

 /*******************************************************************
 * arithmetic_to_boolean: Transforms the given arithmetic share into three
                          boolean shares such that:
                          x = s1 ^ s2 ^ s3
 * *****************************************************************/
 void arithmetic_to_boolean(AShare *x, BShare *s1, BShare *s2, BShare *s3);

/*******************************************************************
 * generate_random_seed: Generates a random seed
 * *****************************************************************/
Seed generate_random_seed();

/*******************************************************************
 * generate_rand_bit_shares: Generates len random bits and their
 *          corresponding arithmetic and binary shares.
 *          This method is used during pre-processing to enable
 *          the binary-to-arithmetic conversion of single bits
 *          required for COUNT().
 * *****************************************************************/
void generate_rand_bit_shares(BShare *rb1, BShare *ra1, 
                              BShare *rb2, BShare *ra2,
                              BShare *rb3, BShare *ra3, int len);

#endif
