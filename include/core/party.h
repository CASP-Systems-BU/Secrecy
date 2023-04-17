#ifndef PARTY_H
#define PARTY_H

#include "mpctypes.h"

/*******************************************************************
 * exchange_rsz_seeds: Initializes the party's random generator for
 *              sharing of 0 and communicates with the other prties.
 *              Each party receives a seed from its predecessor and
 *              send a seed to its successor.
 * 
 *              NOTE: This method must be called before the first 
 *              call to get_next_*()
 * *****************************************************************/
int exchange_rsz_seeds(int succ_rank = -1, int pred_rank = -1);

/*******************************************************************
 * get_next_w: Generates or retrieves the next pair of shares
 *             for the random w value used in arithmetic 
 *             multiplication and equality.
 * *****************************************************************/
WSharePair get_next_w();

/*******************************************************************
 * get_next_r:  Retrieves the next zero-shared arithmetic r value.
 *              
 *              NOTE: exchange_rsz_seeds() must be called before 
 *              the first call to this method.
 * *****************************************************************/
AShare get_next_r();

/*******************************************************************
 * get_next_rb: Retrieves the next zero-shared binary r value.
 * 
 *              NOTE: exchange_rsz_seeds() must be called before 
 *              the first call to this method.
 * *****************************************************************/
BShare get_next_rb();

/*******************************************************************
 * get_next_rb_array: Retrieves len next zero-shared binary r values.
 * 
 *              NOTE: exchange_rsz_seeds() must be called before 
 *              the first call to this method.
 * *****************************************************************/
void get_next_rb_array(BShare *rnum, int len);

/*******************************************************************
 * get_next_array: Retrieves len next zero-shared arithmetic r values.
 * 
 *              NOTE: exchange_rsz_seeds() must be called before 
 *              the first call to this method.
 * *****************************************************************/
void get_next_array(BShare *rnum, int len);

/*******************************************************************
 * get_next_rb_pair_array: Generates an array of random boolean
 *              share pairs.These are used when converting from 
 *              arithmetic to binary.
 * 
 *              NOTE: exchange_rsz_seeds() must be called before 
 *              the first call to this method.
 * *****************************************************************/
void get_next_rb_pair_array(BShare *r1, BShare *r2, int len);

#endif
