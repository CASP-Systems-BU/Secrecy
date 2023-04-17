#ifndef COMM_H
#define COMM_H

#include "mpctypes.h"
#include "assert.h"

/*******************************************************************
 * init:    Initializes MPI and retrieves each party's rank.
 *          This method needs to be called before any other exchange
 *          or computation method.
 * *****************************************************************/
void init(int argc, char** argv);

/*******************************************************************
 * close:   Cleanup and finalize MPI.
 * *****************************************************************/
void close();

/*******************************************************************
 * get_rank:    Returns this party's rank.
 *              init() must have been called before this method.
 * *****************************************************************/
int get_rank();

/*******************************************************************
 * exchange_shares: Binary share exchange.
 *                  Each party sends its share s1 to its predecessor
 *                  and waits to receive the remote share s2
 *                  from its successor.
 *
 *                  If x = x1^x2^x3, after the exchange, party Pi has:
 *                  P1: x1, x2
 *                  P2: x2, x3
 *                  P3: x3, x1
 * *****************************************************************/
BShare exchange_shares(BShare s1);

BShare exchange_shares_async(BShare s1);

/*******************************************************************
 * exchange_shares_u: Similar to exchange_shares() but for unsigned long long
 *                  Each party sends its share s1 to its predecessor
 *                  and waits to receive the remote share s2
 *                  from its successor.
 *
 *                  If x = x1^x2^x3, after the exchange, party Pi has:
 *                  P1: x1, x2
 *                  P2: x2, x3
 *                  P3: x3, x1
 * *****************************************************************/
unsigned long long exchange_shares_u(unsigned long long s1);

/*******************************************************************
 * exchange_bit_shares: Single-bit boolean share exchange.
 *                  Each party sends its share s1 to its predecessor
 *                  and waits to receive the remote share s2
 *                  from its successor.
 *
 *                  If x = x1^x2^x3, after the exchange, party Pi has:
 *                  P1: x1, x2
 *                  P2: x2, x3
 *                  P3: x3, x1
 * *****************************************************************/
BitShare exchange_bit_shares(BitShare s1);

/*******************************************************************
 * exchange_shares_array: Binary share array exchange.
 *                  Each party sends its shares to its predecessor
 *                  and waits to receive the remote shares
 *                  from its successor.
 *
 * *****************************************************************/
void exchange_shares_array(BShare* shares1, BShare* shares2, long length);

/*******************************************************************
 * exchange_shares_array_u: Same as exhcnage_shares_array() but for
 *                          unsigned long long
 *
 * *****************************************************************/
void exchange_shares_array_u(unsigned long long* shares1,
                             unsigned long long* shares2, int length);
/*******************************************************************
 * exchange_a_shares_array: Arithmetic share array exchange.
 *                  Each party sends its shares to its predecessor
 *                  and waits to receive the remote shares
 *                  from its successor.
 *
 * *****************************************************************/
void exchange_a_shares_array(AShare* shares1, AShare* shares2, int length);


/*******************************************************************
 * exchange_bit_shares_array: Single-bit share array exchange.
 *                  Each party sends its shares to its predecessor
 *                  and waits to receive the remote shares
 *                  from its successor.
 *
 * *****************************************************************/
void exchange_bit_shares_array(BitShare* shares1, BitShare* shares2, int length);


/*******************************************************************
 * get_succ: Return this party's successor's MPI rank.
 *           Parties are connected in a ring, i.e. get_succ(P1) = P2
 *           and get_succ(P3) = P1.
 * *****************************************************************/
int get_succ();

/*******************************************************************
 * get_pred: Return this party's predecessor's MPI rank.
 *           Parties are connected in a ring, i.e. get_pred(P1) = P3
 *           and get_pred(P3) = P2.
 * *****************************************************************/
int get_pred();

/*******************************************************************
 * open_b:  Reveal the value of a boolean-shared integer.
 *
 *          P1 serves as the learner who receives the shares from
 *          P2, P3, performs the xor operation, and returns
 *          the result. P2, P3 return their original share.
 * *****************************************************************/
Data open_b(BShare s);

/*******************************************************************
 * open_b_array:  Reveal an array of boolean-shared integers.
 *
 *          P1 serves as the learner who receives the shares from
 *          P2, P3, performs the xor operation, and returns
 *          the result. P2, P3 return their original share.
 *
 *          Updates s in place. It must not be called more than once with the
 *          the same array s.
 *
 * *****************************************************************/
void open_b_array(BShare *s, int len, Data res[]);

void open_byte_array(char *s, int len, char res[]);

/*******************************************************************
 * open_bit_array:  Same as open_b_array but for BitShares
 * *****************************************************************/
void open_bit_array(BitShare *s, int len, bool res[]);

Data open_bit(BitShare s);

/*******************************************************************
 * open_a:  Reveal the value of an arithmetic-shared integer.
 *
 *          P1 serves as the learner who receives the shares from
 *          P2, P3, performs the add operation, and returns
 *          the result. P2, P3 return their original share.
 * *****************************************************************/
Data open_a(AShare s);

/*******************************************************************
 * open_array:  Reveal an array of arithmetic-shared bits.
 *
 *          P1 serves as the learner who receives the shares from
 *          P2, P3, performs the add operation, and returns
 *          the result. P2, P3 return their original share.
 * *****************************************************************/
void open_array(AShare *s, int len, Data res[]);


/*******************************************************************
 * open_mixed_array:  Opens an array that contains both boolean and
 *                    arithmetic shares
 *
 *          a contains the indices of the arithmetic shares
 *          al is the size of a
 *          b contains the indices of the booleabn shares
 *          bl is the size of b
 * *****************************************************************/
//
void open_mixed_array(BShare *s, int rows, int cols, Data res[],
                      unsigned* a, int al, unsigned *b, int bl);

/*******************************************************************
 * reveal_b:  Reveal the value of a boolean-shared integer and
 *          share it with every party.
 *
 *          P1 serves receives the shares from P2, P3,
 *          performs the xor operation, and sends back
 *          the result. All P1, P2, P3 learn (return) the value.
 * *****************************************************************/
Data reveal_b(BShare s);

/*******************************************************************
 * reveal_b_array:  Array-based implementation of reveal_b().
 *                  The result is stored in s.
 * *****************************************************************/
void reveal_b_array(BShare *s, int len);
void reveal_b_array_async(BShare *s, int len);

#endif
