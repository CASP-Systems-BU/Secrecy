#ifndef MPCTYPES_H
#define MPCTYPES_H

#include <limits.h>
#include <stdbool.h>

/**
 * Original data
 **/
typedef long long Data;

/**
 * Plaintext table
 **/
typedef Data** DataTable;

/**
 * Seed for random number generator
 **/
typedef long long Seed;

/**
 * Arithmetic Share
 **/
typedef long long AShare;

/**
 * Boolean Share
 **/
typedef long long BShare;

typedef long long Share;
typedef Share *Shares;
typedef AShare *AShares;
typedef BShare *BShares;

/**
 * Boolean single-bit Share
 **/
typedef bool BitShare;

/**
 * represents an share table of an input Table
 * that was produced by a data owner with ownerId
 * and was assigned to the computation party with partyId.
**/
typedef struct {
    int ownerId;
    int partyId;
    int numRows;
    int numCols;
    int relationId;
    Share **content;
} ShareTable, AShareTable, BShareTable;

/**
 * Table represents a relational table provided by the data owner with  ownerId.
**/
typedef struct {
    int ownerId;
    int numRows;
    int numCols;
    Data **content;
} Table;

/**
 * PartyMessage represents a data message exchanges between data owners.
**/
typedef struct {
    int relationId;
    int rowNumber;
    // Payload p; //TODO: what's the type of Payload?
} PartyMessage;


/**
 * A pair of shares for a random w that is used in multiplications and arithmetic equality.
 **/
typedef struct {
    AShare first;
    AShare second;
} WSharePair;

#endif
