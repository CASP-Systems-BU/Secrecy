#ifndef SECRECY_PROJECT_OPERATOR_H
#define SECRECY_PROJECT_OPERATOR_H

// project is used to reduce the number of used columns
// - Instant projection in plan: to reduce memory and done instantly
// - Delay projection in plan: to reduce memory in future steps and delayed to
//  other table cardinality changer (i.e join or another instant project)

/**
 * The ProjectOperator reduces the number of columns in the input data.
 */
class ProjectOperator : public Operator{

    virtual bool internalExecute() {
        output = inputs[0];
        return true;
    }

public:


    /**
     * The ProjectOperator reduces the number of columns in the input data.
     * @param input_1 The previous operator which this operator will process its output.
     * @param expressions Expressions to specify the columns to be kept.
     * @param batchSize Batch size this operator will use to process input data.
     */
    ProjectOperator(std::shared_ptr<Operator> input_1,
                   const std::vector<std::shared_ptr<Expression>> &expressions,
                   const int &batchSize)
    : Operator(OperatorType::Projection, batchSize) {
        this->operators.push_back(input_1);
        this->expressions = expressions;
    }


    /**
     * The ProjectOperator reduces the number of columns in the input data.
     * Note: This constructor is used initialize the expressions and batchSize.
     *  However, the input operator or table has to be modified before execution.
     * @param expressions Expressions to specify the columns to be kept.
     * @param batchSize Batch size this operator will use to process input data.
     */
    ProjectOperator(const std::vector<std::shared_ptr<Expression>> &expressions,
                   const int &batchSize)
            : Operator(OperatorType::Projection, batchSize) {
        this->expressions = expressions;
    }

};

#endif //SECRECY_PROJECT_OPERATOR_H