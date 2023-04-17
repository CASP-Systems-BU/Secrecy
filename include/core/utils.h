#ifndef UTILS_H
#define UTILS_H

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <sodium.h>
#include "mpi.h"

#include "mpctypes.h"
#include "sharing.h"

/*******************************************************************
* init_tables: Initializes shares tables and distributes random shares
*              across parties
* *****************************************************************/
void init_tables(BShareTable *t1, BShareTable *t2);

/*******************************************************************
 * allocate_int_shares_table: Allocates 2D array of given size to store shares
 *
 * *****************************************************************/
void allocate_int_shares_table(AShareTable *table);

/*******************************************************************
 * allocate_bool_shares_table: Allocates 2D array of given size to
 *                             store boolean shares
 * *****************************************************************/
void allocate_bool_shares_table(BShareTable *table);

/*******************************************************************
 * allocate_a_shares_table: Allocates 2D array of given size to
 *                             store arithmetic shares
 * *****************************************************************/
void allocate_a_shares_table(AShareTable *table);

/*******************************************************************
 * populate_shares: Populates shares table
 *
 * *****************************************************************/
void populate_shares_table(AShareTable* party_shares, AShare* all_shares,
                           int columns, int* share_indices);

/*******************************************************************
* print_shares_table: Prints table of shares to standard output
*
* *****************************************************************/
void print_shares_table(AShareTable* party_shares);

/*******************************************************************
* print_bool_shares_table: Prints table of shares to standard output
*
* *****************************************************************/
void print_bool_shares_table(BShareTable* party_shares);


/*******************************************************************
* generate_random_table: Generates a random table of a given size
* *****************************************************************/
void generate_random_table(Table *table, int rows, int columns);

// TODO: inline (get_bit - get_bit_u - get_bit_u8 - unset_lsbs - to_bshare)

/*******************************************************************
* get_bit: Returns the i-th bit the given boolean share
*          as a boolean share
* *****************************************************************/
BShare get_bit(const BShare s, int i);

/*******************************************************************
* get_bit_u: Returns the i-th bit the given unsigned long long
* *****************************************************************/
BShare get_bit_u(unsigned long long s, int i);

BShare get_bit_u8(char s, char i);

/*******************************************************************
* set_bit: Sets the i-th bit
*          and returns a new boolean share
* *****************************************************************/
BShare set_bit(const BShare s, int i);

/*******************************************************************
* unset_bits: Unsets the i LSBs of the given boolean share
*             and returns a new boolean share
* *****************************************************************/
BShare unset_lsbs(const BShare s, int i);

/*******************************************************************
* unset_bits: Unsets the i-th bit of the given boolean share
*             and returns a new boolean share
* *****************************************************************/
BShare unset_bit(const BShare s, int i);

/*******************************************************************
* to_bshare: Returns a boolean share equivlent to the given single-bit boolean
*            share, i.e., if s is set, it returns 0xF....F, else returns 0
* *****************************************************************/
BShare to_bshare(const BitShare s);

/*******************************************************************
* print_binary: Prints a BShare in binary format
* *****************************************************************/
void print_binary(const BShare s);

Data** allocate_table(long rows, long columns);

BShare** allocate_2D_table(long numRows, long numCols);

Data** allocate_2D_data_table(long numRows, long numCols);

BitShare** allocate_2D_bit_table(long numRows, long numCols);

char** allocate_2D_byte_array(long numRows, long numCols);

int** allocate_int_2D_table(long numRows, long numCols);

void generate_and_share_random_data(int rank, BShare *r1s1, BShare *r1s2, long ROWS);

#endif
