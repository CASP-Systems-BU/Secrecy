#ifndef SECRECY_SELECT_OPERATOR_H
#define SECRECY_SELECT_OPERATOR_H

#include "operator.h"


/**
 * The SelectOperator is used to apply a filter to its input by updating
 * the padded Boolean Selection Bit "[SEL]".
 * Note: the padded Arithmetic Selection Bit "SEL" will be corrupted.
 */
class SelectOperator : public Operator {

    // Executes the filter
    virtual bool internalExecute() {
        output = inputs[0];

        // All predicates are binary
        std::shared_ptr<Expression> left = filterExpr->expressions[0]; // left subexpression
        std::shared_ptr<Expression> right = filterExpr->expressions[1]; // right subexpression

        // Get selection predicate type
        switch (filterExpr->expressionType) {
            case Expression::Equal: {
                // Right subexpression is constant
                if (right->expressionType == Expression::LiteralInt) {
                    // get left column index
                    int left_pos = output[left->name].getIndex();
                    // Convention: indices of att-c and c-att columns are next
                    int left_ind = (left_pos+1) * CR_P.PARTY_REPLICATION;
                    int right_ind = (left_pos+2) * CR_P.PARTY_REPLICATION;
                    Predicate_B COND = {EQ, NULL, NULL, left_ind, right_ind};

                    // Create a column for the result
                    DerivedColumn<BShare, BShareTable> sel_c(output.getSize());
                    auto sel_column = sel_c.getShares();
                    auto sel_column_ptr = sel_column.get();

                    select_b(output.shareTable, COND, &sel_column_ptr[0][0]);

                    // Exchange communication
                    exchangeArrayWithParties<BShare>(sel_column_ptr[0], sel_column_ptr[1], output.getSize(), 1,
                                            PRIMITIVES_SHARE_TAG, 1);

                    // append the selection result column to the table
                    output["[SEL]"] = sel_c;
                    output.bSelectionUsed = true;
                }
                else {
                    std::cout << "EXPRESSION LEFT " << left->expressionType  << " RIGHT " <<  right->toString() << std::endl;
                    assert((right->expressionType == Expression::AttributeRef ||
                            right->expressionType == Expression::DerivedAttributeRef) &&
                           (left->expressionType == Expression::AttributeRef ||
                                   left->expressionType == Expression::DerivedAttributeRef));

                    if(output.bSelectionUsed){
                        output["[SEL]"] = output["[SEL]"] & (output[right] == output[left]);
                    }else{
                        output["[SEL]"] = (output[right] == output[left]);
                        output.bSelectionUsed = true;
                    }
                    output.aSelectionUsed = false;

                    break;
                }
                break;
            }
            case  Expression::LessThanEqual:
            {
                if (right->expressionType == Expression::LiteralInt) {
                    leaderLogging("[ERROR] Unsupported filter type.\n", true);
                }else {
                    assert(right->expressionType == Expression::AttributeRef && left->expressionType == Expression::AttributeRef);

                    if (output.bSelectionUsed) {
                        output["[SEL]"] = output["[SEL]"] & (output[right] <= output[left]);
                    } else {
                        output["[SEL]"] = (output[right] <= output[left]);
                        output.bSelectionUsed = true;
                    }
                    output.aSelectionUsed = false;
                }

                break;
            }
            default:
                leaderLogging("[ERROR] Unknown filter predicate: " + (filterExpr->name) + "\n", true);
                break;
        }
        return true;
    }

public:
    /***
     * The filter expression
     */
    std::shared_ptr<Expression> filterExpr;

    /***
     * The SelectOperator is used to apply a filter to its input by updating
     * the padded Boolean Selection Bit "[SEL]".
     * Note: the padded Arithmetic Selection Bit "SEL" will be corrupted.
     * @param input_1 - The previous operator which this operator will process its output.
     * @param expressions - Expressions to specify a set of filtering conditions.
     * @param batchSize - Batch size this operator will use to process input data.
     */
    SelectOperator(std::shared_ptr<Operator> input_1, const std::vector<std::shared_ptr<Expression>> &expressions,
                   const int &batchSize)
        : Operator(OperatorType::Selection, batchSize)
    {
        this->operators.push_back(input_1);
        this->expressions = expressions;
    }


    /**
     * The SelectOperator is used to apply a filter to its input by updating
     * the padded Boolean Selection Bit "[SEL]".
     * Note: the padded Arithmetic Selection Bit "SEL" will be corrupted.
     * @param input_1 - The previous operator which this operator will process its output.
     * @param filterExpr - Expressions to specify the filtering condition.
     * @param batchSize - Batch size this operator will use to process input data.
     */
    SelectOperator(std::shared_ptr<Operator> input_1, const std::shared_ptr<Expression> &filterExpr,
                   const int &batchSize)
        : Operator(OperatorType::Selection, batchSize)
    {
        this->operators.push_back(input_1);
        this->filterExpr = filterExpr;
    }

    /***
     * The SelectOperator is used to apply a filter to its input by updating
     * the padded Boolean Selection Bit "[SEL]".
     * Note: the padded Arithmetic Selection Bit "SEL" will be corrupted.
     * Note: This constructor is used initialize the expressions and batchSize.
     *  However, the input operator or table has to be modified before execution.
     * @param expressions - Expressions to specify a set of filtering conditions.
     * @param batchSize - Batch size this operator will use to process input data.
     */
    SelectOperator(const std::vector<std::shared_ptr<Expression>> &expressions, const int &batchSize)
        : Operator(OperatorType::Selection, batchSize)
    {
        this->expressions = expressions;
    }

};


#endif //SECRECY_SELECT_OPERATOR_H