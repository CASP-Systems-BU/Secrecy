#ifndef SECRECY_TABLE_H
#define SECRECY_TABLE_H

#include <iostream>
#include "constants.h"
#include "../core/mpctypes.h"
#include "../core/sharing.h"

#include <memory>

typedef std::shared_ptr<BShares> RowData;
typedef std::shared_ptr<BShares> ColumnData;

// Assisting function for memory allocation
template<typename T>
T **allocateContent(int rows, int cols, int replicas = CR_P.PARTY_REPLICATION, bool zeros = false) {
    long long rowsPart = rows * sizeof(T *);
    long long colsPart = sizeof(T) * rows * cols * replicas;
    long long length = rowsPart + colsPart;

    T **content;

    if (zeros) {
        length = length / (sizeof(char));
        content = (T **) calloc(length, sizeof(char));
    } else {
        content = (T **) malloc(length);
    }
    assert(content != nullptr);

    T *ptr = (T *) (content + rows);

    for (int i = 0; i < rows; i++) {
        content[i] = (ptr + cols * replicas * i);
    }

    return content;
}

template<typename T, typename T2>
void allocateContent(T2 &table, int replicas = CR_P.PARTY_REPLICATION, bool zeros = true) {
    table.content = allocateContent<T>(table.numRows, table.numCols, replicas, zeros);
}


template<typename T>
std::shared_ptr<T *>
allocateColumnWiseContent(int rows, int cols, int replicas = CR_P.PARTY_REPLICATION, bool zeros = false) {
    int colsPart = cols * replicas * sizeof(T *);
    int rowsPart = cols * replicas * rows * sizeof(T);
    int length = colsPart + rowsPart;

    T **content_ptr;

    if (zeros) {
        length = length / (sizeof(T **)) + (length % sizeof(T **) > 0 ? 1 : 0);
        content_ptr = (T **) calloc(length, sizeof(char));
    } else {
        content_ptr = (T **) malloc(length);
    }
    assert(content_ptr != nullptr);

    T *ptr = (T *) (content_ptr + cols * replicas);

    for (int i = 0; i < cols * replicas; i++) {
        content_ptr[i] = (ptr + rows * i);
    }

    auto content = std::shared_ptr<T *>(content_ptr);

    return content;
}

// Share Table
// - Getting inputs.
// - From another share table.
// - from a Data table.
// Parameters: rows, cols, [content, shareTable, dataTable],{ownerID, partyID, relationID}
// num of cols used is the number of attributes
template<typename T, typename T2>
T2 createNewShareTable(int rows, int cols, int partyID = -1,
                       int relationID = -1, int ownerID = -1) {

    // create the share table Object
    T2 shareTable = T2{ownerID, partyID, rows, cols, relationID};
    allocateContent<T, T2>(shareTable);

    return shareTable;
}


template<typename T, typename T2>
T2 newShareTable(int rows, int cols, T **content, int partyID = -1, int relationID = -1, int ownerID = -1) {

    T2 shareTable = T2{ownerID, partyID, rows, cols, relationID};
    shareTable.content = content;

    return shareTable;
}

template<typename T>
void freeShareTable(T *table) {
    if (table != nullptr) {
        if (table->content != nullptr) {
            free(table->content);
        }
        free(table);
    }
}

template<typename T>
void freeShareTable(T table) {
    if (table.content != nullptr) {
        free(table.content);
    }
}


template<typename T, typename T2>
void addToTable(T2 &bShareTable, std::shared_ptr<T *> &shares, const int &col_1, const int &col_2) {
    auto shares_ptr = shares.get();
    // Fill the output table
    for (int i = 0; i < bShareTable.numRows; i++) {
        bShareTable.content[i][col_1] = shares_ptr[0][i];
        bShareTable.content[i][col_2] = shares_ptr[1][i];
    }
}

template<typename T, typename T2>
void addToTable(T2 &bShareTable, std::shared_ptr<T *> &shares, const int &col_1) {
    int col__1 = col_1 * R_3PC.PARTY_REPLICATION;
    addToTable<T, T2>(bShareTable, shares, col__1, col__1 + 1);
}


template<typename T, typename T2>
void addToTable(T2 &bShareTable, T **shares, const int &col_1, const int &col_2) {
    // Fill the output table
    for (int i = 0; i < bShareTable.numRows; i++) {
        bShareTable.content[i][col_1] = shares[0][i];
        bShareTable.content[i][col_2] = shares[1][i];
    }
}

template<typename T, typename T2>
void addToTable(T2 &bShareTable, T **shares, const int &col_1) {
    int col__1 = col_1 * R_3PC.PARTY_REPLICATION;
    addToTable<T, T2>(bShareTable, shares, col__1, col__1 + 1);
}

template<typename T, typename T2>
std::shared_ptr<T *> getFromTable(const T2 &bShareTable, int col_1) {
    assert(col_1 < bShareTable.numCols);

    col_1 *= R_3PC.PARTY_REPLICATION;
    int col_2 = col_1 + 1;
    auto shares = allocateColumnWiseContent<T>(bShareTable.numRows, 1);
    auto shares_ptr = shares.get();

    for (int i = 0; i < bShareTable.numRows; ++i) {
        shares_ptr[0][i] = bShareTable.content[i][col_1];
        shares_ptr[1][i] = bShareTable.content[i][col_2];
    }

    return shares;
}

// Validation
template<typename T>
void validateTableInput(const T &inputOne, const T &inputTwo,
                        const T &output, const int &colOne,
                        const int &colTwo, const int &colOut) {

    assert(inputOne.numRows == inputTwo.numRows);
    assert(inputOne.numRows == output.numRows);
    assert(inputOne.numCols > colOne);
    assert(inputTwo.numCols > colTwo);
    assert(output.numCols > colOut);
}

template<typename T>
void validateTableInput(const T &inputOne, const T &output,
                        const int &colOne, const int &colOut) {

    assert(inputOne.numRows == output.numRows);
    assert(inputOne.numCols > colOne);
    assert(output.numCols > colOut);
}

#endif //SECRECY_TABLE_H
