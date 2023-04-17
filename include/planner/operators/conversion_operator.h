#ifndef SECRECY_CONVERSION_OPERATOR_H
#define SECRECY_CONVERSION_OPERATOR_H

#include "operator.h"

/**
 * The ConversionOperator converts From a type of secret sharing to another.
 * Conversion type is deduced from input and output types.
 * For example, using column "A" as input and column "[B]" as an output will make automatic
 * conversion from arithmetic secret sharing to boolean secret sharing.
 * Currently, supporting only "A" to "[B]" conversion.
 */
class ConversionOperator : public Operator {
    std::vector<std::shared_ptr<Expression>> dest_expressions;

    virtual bool internalExecute() {
        output = inputs[0];

        // Currently, supporting only "A" to "[B]" conversion.
        for(int i =0; i < this->expressions.size(); ++i){
            output[this->dest_expressions[i]] = convertToBColumn<BShare, BShareTable>(output[this->expressions[i]]);
        }

        return true;
    }


public:

    /**
     * The ConversionOperator converts From a type of secret sharing to another.
     * Conversion type is deduced from input and output types.
     * For example, using column "A" as input and column "[B]" as an output will make automatic
     * conversion from arithmetic secret sharing to boolean secret sharing.
     * Currently, supporting only "A" to "[B]" conversion.
     * @param input_1 The previous operator which this operator will process its output.
     * @param expressions Expressions to specify the input columns.
     * @param dest_expressions Expressions to specify the output columns.
     * @param batchSize Batch size this operator will use to process input data.
     */
    ConversionOperator(std::shared_ptr<Operator> input_1,
                       const std::vector<std::shared_ptr<Expression>> &expressions,
                       const std::vector<std::shared_ptr<Expression>> &dest_expressions,
                       const int batchSize = 8)
            : Operator(OperatorType::Conversion, batchSize) {
        this->operators.push_back(input_1);
        this->expressions = expressions;
        this->dest_expressions = dest_expressions;
    }

};

#endif //SECRECY_CONVERSION_OPERATOR_H
