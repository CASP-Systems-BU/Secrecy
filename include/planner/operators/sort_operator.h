#ifndef SECRECY_SORT_OPERATOR_H
#define SECRECY_SORT_OPERATOR_H

#include "operator.h"

/**
 * The SortOperator sorts a given table using the data in a few columns in ascending or descending order.
 * Note: Arithmetic shares in this table will be corrupted.
 */
class SortOperator : public Operator {
    std::vector<bool> asc;

    virtual bool internalExecute() {
        output = inputs[0];
        unsigned num_keys = expressions.size();
        bool asc_[num_keys];
        unsigned key_pos[num_keys];

        // get key positions
        for (int i=0; i<num_keys; i++) {
            key_pos[i] = output[expressions[i]->name].getIndex();
        }

        for (int i = 0; i < num_keys; ++i) {
            asc_[i] = asc[i];
            key_pos[i] = key_pos[i] * CR_P.PARTY_REPLICATION;
        }

        if (batchSize == output.getSize()) {
            batchSize /= 2;
        }
        output.shareTable.numCols *= CR_P.PARTY_REPLICATION;
        bitonic_sort_batch(&output.shareTable, key_pos, num_keys, asc_, batchSize);
        output.shareTable.numCols /= CR_P.PARTY_REPLICATION;

        return true;
    }


public:

    /**
     * The SortOperator sorts a given table using the data in a few columns in ascending or descending order.
     * Note: Arithmetic shares in this table will be corrupted.
     * @param input_1 The previous operator which this operator will process its output.
     * @param expressions Expressions to specify the columns used for sorting in order.
     * @param asc Booleans to specify whether each used column is sorted in ascending order.
     * @param batchSize Batch size this operator will use to process input data
     */
    SortOperator(std::shared_ptr<Operator> input_1, const std::vector<std::shared_ptr<Expression>> &expressions,
                 const std::vector<bool> &asc, const int &batchSize)
            : Operator(OperatorType::Sort, batchSize) {
        this->operators.push_back(input_1);
        this->expressions = expressions;
        this->asc = asc;
    }

};


#endif //SECRECY_SORT_OPERATOR_H