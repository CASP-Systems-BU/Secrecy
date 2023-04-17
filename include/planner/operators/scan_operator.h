#ifndef SECRECY_SCAN_OPERATOR_H
#define SECRECY_SCAN_OPERATOR_H

#include "operator.h"

/**
 * The ScanOperator is designed to read data from a relation in batches.
 */
class ScanOperator : public Operator {

    virtual bool internalExecute() {
        output = inputs[0];  // Relation is already in memory
        for (int j=0; j<output.getCols().size(); j++) {
            if (output.getCols()[j] == B_SEL) {
                output.bSelectionUsed = true;
                output.initBColumnOnes(j);  // Initialize with 1's (1^1^1=1)
                break;
            }
        }
        return true;
    }

    virtual std::string getInternalString() {
        return inputs[0].getName();
    }

public:
    /**
     * The ScanOperator is designed to read data from a relation in batches.
     * @param relation The input relation that has the data to be processed.
     * @param batchSize Batch size this operator will use to process input relation data.
     */
    ScanOperator(const std::shared_ptr<Relation> &relation, int batchSize = 8) :
        Operator(OperatorType::Scan, batchSize)
    {
        inputs.push_back(*relation);
    }

    /**
     * The ScanOperator is designed to read data from a relation in batches.
     * Note: This constructor is used initialize the batchSize.
     *  However, the input operator or table has to be modified before execution.
     * @param relation The input relation that has the data to be processed.
     * @param batchSize Batch size this operator will use to process input relation data.
     */
    ScanOperator(int batchSize = 8) : Operator(OperatorType::Scan, batchSize) {}
};


#endif //SECRECY_SCAN_OPERATOR_H