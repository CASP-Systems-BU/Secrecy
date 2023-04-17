#ifndef SECRECY_COMMUNICATION_H
#define SECRECY_COMMUNICATION_H

#include "constants.h"
#include "../core/comm.h"
#include "table.h"

// Some predefined tags
#define TEST_SHARE_TAG 14
#define PRIMITIVES_SHARE_TAG 21
#define Distinct_SHARE_TAG 7

// Send array to party
#include <type_traits>

// TODO: create another inner template functions in .cpp
// - Can the header and implementation of the template function be seperated?
// - To remove dependency of the MPI on upper layers.
#include "mpi.h"

// TODO: refactor to have this format Communication_Round(T** data, int rows)

template<typename T>
void sendArrayToParty(T *data, int rows, int cols, int partyID, int shareTag,
                      int replicas = CR_P.PARTY_REPLICATION) {
    int length = rows * cols * replicas;
    if ((std::is_same<T, long long>::value)) {
        MPI_Send(data, length, MPI_LONG_LONG, partyID, shareTag, MPI_COMM_WORLD);
    }
}

template<typename T>
void getArrayFromParty(T *data, int rows, int cols, int partyID, int shareTag,
                       int replicas = CR_P.PARTY_REPLICATION) {
    int length = rows * cols * replicas;
    if ((std::is_same<T, long long>::value)) {
        MPI_Recv(data, length, MPI_LONG_LONG, partyID, shareTag, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }
}

template<typename T>
void exchangeArrayWithParties(T *dataOut, T *dataIn,
                              int rows, int cols,
                              int outPartyID, int inPartyID,
                              int shareTag, int replicas = CR_P.PARTY_REPLICATION) {
    sendArrayToParty<T>(dataOut, rows, cols, outPartyID, shareTag);
    getArrayFromParty<T>(dataIn, rows, cols, inPartyID, shareTag);
}

template<typename T>
void exchangeArrayWithParties(T *dataOut, T *dataIn,
                              int rows, int cols,
                              int shareTag, int replicas = CR_P.PARTY_REPLICATION) {
    if ((std::is_same<T, long long>::value)) {
        exchange_shares_array(dataOut, dataIn, rows * cols * replicas);
    } else {
        int fromParty = get_succ();
        int toParty = get_pred();

        sendArrayToParty<T>(dataOut, rows, cols, toParty, shareTag, replicas);
        getArrayFromParty<T>(dataIn, rows, cols, fromParty, shareTag, replicas);
    }
}


// Table assisting function
template<typename T, typename T2>
void exchangeAndAddToTable(T2 &bShareTable, T **shares, const int &col_1, const bool &freeShares = true) {
    // Communication round
    int fromParty = get_succ();
    int toParty = get_pred();
    exchangeArrayWithParties<T>(shares[0], shares[1], bShareTable.numRows, 1, toParty, fromParty,
                                PRIMITIVES_SHARE_TAG, 1);

    // Fill the output table
    addToTable<T, T2>(bShareTable, shares, col_1);

    // Free temp arrays;
    if (freeShares) {
        free(shares);
    }
}

template<typename T, typename T2>
void exchangeAndAddToTable(T2 &bShareTable, std::shared_ptr<T *> &shares, const int &col_1, const bool &freeShares = true) {
    auto shares_ptr = shares.get();
    // Communication round
    int fromParty = get_succ();
    int toParty = get_pred();
    exchangeArrayWithParties<T>(shares_ptr[0], shares_ptr[1], bShareTable.numRows, 1, toParty, fromParty,
                                PRIMITIVES_SHARE_TAG, 1);

    // Fill the output table
    addToTable<T, T2>(bShareTable, shares, col_1);
}


#endif //SECRECY_COMMUNICATION_H
