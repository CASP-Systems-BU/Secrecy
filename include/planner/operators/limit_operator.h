#ifndef SECRECY_LIMIT_OPERATOR_H
#define SECRECY_LIMIT_OPERATOR_H

#include "operator.h"

/**
 * The LimitOperator reduces the cardinality of its input to a certain size.
 * Note: If the table's size is already smaller than the limitSize, no changes
 * happen to input table.
 */
class LimitOperator : public Operator {

    // TODO: implement this operator.
    virtual bool internalExecute() {
        output = inputs[0];
        return true;
    }

public:

    /**
     * The LimitOperator reduces the cardinality of its input to a certain size.
     * Note: If the table's size is already smaller than the limitSize, no changes
     * happen to input table.
     * @param input_1 The previous operator which this operator will process its output.
     * @param limitSize The max size the input data will trimmed to.
     * @param batchSize Batch size this operator will use to process input data.
     */
    LimitOperator(std::shared_ptr<Operator> input_1,
                  const int &limitSize,
                  const int &batchSize)
    : Operator(OperatorType::Limit, batchSize) {
        this->operators.push_back(input_1);
        std::shared_ptr<Expression> limit(new ConstantExpr((int64_t) limitSize));
        this->expressions.push_back(limit);
    }
};


#endif //SECRECY_LIMIT_OPERATOR_H
