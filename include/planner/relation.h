#ifndef SECRECY_RELATION_H
#define SECRECY_RELATION_H

#include "../common/table.h"
#include "table/column.h"
#include "table/derived_column.h"
#include "table/a_column.h"
#include "table/b_column.h"
#include "table/column_conversion.h"
#include "operators/expression.h"

#include <string>
#include <vector>
#include <map>
#include <iostream>

// TODO: Move to utils
/**
 * Prints a string from the computing party with ID=0.
 * @param str - The string to print.
 * @param error - Whether to exit the program returning an error code.
 */
inline void leaderLogging(std::string str, bool error=false) {
    if (get_rank() == 0) {
        std::cout << str;
    }
    if (error)
        exit(-1);
}

/**
 * Prints a string followed by an array elements separated by tab
 * from the computing party with ID=0.
 * @param str - The string to print.
 * @param array - The array to print its elements.
 * @param a_size - The size of the elements in the array.
 * @param error - Whether to exit the program returning an error code.
 */
inline void leaderLogging(std::string str, int *array, int a_size, bool error=false) {
    if (get_rank() == 0) {
        for (int i = 0; i < a_size; i++)
            str += std::to_string(array[i]) + "  ";
        str += '\n';
        std::cout << str;
    }
    if (error)
        exit(-1);
}

typedef std::map<std::string, Column *> Schema;

class Relation {
    // Relation's unique id
    int id;
    // Schema
    Schema schema;
    // Metadata
    std::string name;
    std::string versionName;
    bool materialized;

    Column &getColumnA(const std::string &columnName);
    Column &getColumnB(const std::string &columnName);

public:
    // Data
    ShareTable shareTable;

    // Relation alias
    std::string alias;

    // Col names
    std::vector<std::string> colNames;
    // Boolean shares
    std::vector<bool> bShares;

    bool bSelectionUsed = false;
    bool aSelectionUsed = false;

    Relation();

    Relation(const std::string &name, const std::vector<std::string> &cols, const int &rows);

    Relation(const std::string &name, const std::vector<std::string> &cols, const int &rows, std::vector<bool> &bShares);

    ~Relation();

    /**
     * Allocate data memory if not already allocated.
     * Assigns each column in the schema to its
     * corresponding memory location.
     */
    void allocateTable();

    /**
     * Frees the assigned memory for data if already allocated.
     */
    void freeRelation();

    /**
     * Adds a new column to the schema and point it to
     * the first corresponding not assigned column in the table.
     * Note: adding "[]" to column name means it's Boolean.
     * Note: does not increase already allocated memory.
     * @param cols - A vector containing the columns names.
     */
    void addColumns(const std::vector<std::string> &cols);

    /**
     * Adds a new column to the schema and point it to
     * the first corresponding not assigned column in the table.
     * Note: adding "[]" to column name means it's Boolean.
     * Note: does not increase already allocated memory.
     * @param col the column name.
     */
    void addColumns(const std::string &col);

    /**
     * The function allocates the needed data memory for the relation
     * and overwrite it to match the given ShareTable.
     * @param other - The allocated ShareTable to copy from.
     */
    void insertTable(const ShareTable &other);

    /**
     * The access operator take s a column name as a string and gives
     * back reference to the corresponding Column that has the secret shares.
     * @param columnName - The name of the requested column.
     * @return - The reference to the column secret shares.
     */
    Column &operator[](const std::string &columnName);

    /**
     * The access operator take s an expression and gives
     * back reference to the corresponding Column that
     * has the secret shares.
     * @return
     */

    /**
     * The access operator take s an expression and gives
     * back reference to the corresponding Column that
     * has the secret shares.
     * @param columnExpression - An expression to specify the requested column.
     * @return - The reference to the column secret shares.
     */
    Column &operator[](const std::shared_ptr<Expression> columnExpression);

    /**
     * Get the relation name as a string.
     * @return - String having relation name.
     */
    std::string getName() const;

    /**
     * Given a column name, decide whether the name should
     * refer to a Boolean Secret Shares Column.
     * @param columnName The column name to test.
     * @return returns
     */
    bool isBColumnName(const std::string &columnName) const;

    /**
     * Returns the number of rows in the relation. All columns
     * have the same number of rows.
     * @return - The number of rows.
     */
    int getSize() const;


    /**
     * Returns a vector containing the names of
     * the columns with Arithmetic secret shares.
     * @return A vector of strings containing the names.
     */
    std::vector<std::string> getACols() const;

    /**
     * Returns a vector containing the names of
     * the columns with Boolean Secret Shares.
     * @return A vector of strings containing the names.
     */
    std::vector<std::string> getBCols() const;

    /**
     * Returns the names of all the columns in the relation.
     * @param omitTableName - Whether to add the table name at the begining.
     * @return A vector of strings containing the names.
     */
    std::vector<std::string> getCols(bool omitTableName=false) const;

    /**
     * Populate the column at the given index with secret shares
     * that make the data stored in the Arithmetic Secret Shares
     * column equal to one.
     * @param index - The column to be modified index.
     */
    void initAColumnOnes(int index);

    /**
     * Populate the column at the given index with secret shares
     * that make the data stored in the Arithmetic Secret Shares
     * column equal to zero.
     * @param index - The column to be modified index.
     */
    void initAColumnZeros(int index);

    /**
     * Populate the column at the given index with secret shares
     * that make the data stored in the Boolean Secret Shares
     * column equal to one.
     * @param index - The column to be modified index.
     */
    void initBColumnOnes(int index);


    /**
     * Returns a vector with the indices of the columns specified.
     * @param colExpressions - A vector of Expressions to specify the columns.
     * @return A pair of an array pointer containing the indices and the number of indices.
     */
    std::pair<unsigned int *, int> getIndices(const std::vector<std::shared_ptr<Expression>> &colExpressions);
    std::shared_ptr<Data *> openRelation() const;


    /**
     * Returns the number of rows in the relation. All columns
     * have the same number of rows.
     * @return - The number of rows.
     */
    unsigned cardinality() { return shareTable.numRows;};

    /**
     * Returns the number of columns in the relation.
     * @return - The number of columns.
     */
    unsigned colsNum() const {return schema.size(); };

    friend class PlanGen;
};

#endif //SECRECY_RELATION_H
