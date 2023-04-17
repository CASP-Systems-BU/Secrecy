#ifndef SECRECY_SEMI_JOIN_OPERATOR_H
#define SECRECY_SEMI_JOIN_OPERATOR_H

#include "operator.h"
#include "../../core/relational.h"


/**
 * The SemiJoinOperator updates the Boolean Selection Bit "[SEL]" for the left table
 *  based on an expression evaluation involving another table. The operator also applies
 *  partial aggregations if any specified.
 */
class SemiJoinOperator : public Operator {

    virtual bool internalExecute() {
        auto r1 = inputs[0];
        auto r2 = inputs[1];
        output = inputs[0];

        Column &c1 = r1[expressions[0]->expressions[0]->toString()];
        Column &c2 = r2[expressions[0]->expressions[1]->toString()];

        if (partialAgg.size() > 0) {
            int rows2 = r2.getSize();
            auto agg = partialAgg[0]; // function expression
            if (agg->name == "COUNT" && agg->expressions[0]->hasType(Expression::Star)) {

                // add right_att "SEL" to right relation
                std::vector<std::string> cols = r2.getCols();
                cols.push_back("SEL");
                auto new_r2 = Relation(r2.getName(), cols, rows2);
                new_r2.allocateTable();
                // copy r2 columns
                int base_col = CR_P.PARTY_REPLICATION * (cols.size()-1);
                for (int i = 0; i < rows2; ++i) {
                    for (int j = 0; j < base_col; ++j) {
                        new_r2.shareTable.content[i][j] = r2.shareTable.content[i][j];
                    }
                }
                int right_att_index = new_r2["SEL"].getIndex();

                // check if right input has [SEL] column
                // and if yes, convert it to arithmetic
                if (r2.bSelectionUsed) {
                    //convert it to arithmetic
                    new_r2["SEL"] = convertToAColumn<AShare, AShareTable>(r2["[SEL]"]);
                }
                // if no, add one with 1s
                else {
                    new_r2.initAColumnOnes(right_att_index);
                }

                // add arithmetic sum column to left input
                std::vector<std::string> cols_left = r1.getCols();
                cols_left.push_back("SEL");
                auto new_r1 = Relation(r1.getName(), cols_left, r1.getSize());
                new_r1.allocateTable();

                // copy r1 columns
                base_col = CR_P.PARTY_REPLICATION * (cols_left.size()-1);
                for (int i = 0; i < r1.getSize(); ++i) {
                    for (int j = 0; j < base_col; ++j) {
                        new_r1.shareTable.content[i][j] = r1.shareTable.content[i][j];
                    }
                }
                unsigned left_sum_index = new_r1["SEL"].getIndex();
                new_r1.initAColumnZeros(left_sum_index);

                int num_batches = r1.getSize() / batchSize;
                AShare *r_a = new AShare[rows2*batchSize];
                BShare *r_b = new BShare[rows2*batchSize];
                get_next_rb_array(r_b, rows2*batchSize);
                get_next_array(r_a, rows2*batchSize);

                for (int i=0; i<num_batches; i+=batchSize) {
                    group_by_join_first(&new_r1.shareTable, &new_r2.shareTable, i, i+batchSize, 
                                c1.getIndex() * CR_P.PARTY_REPLICATION,
                                c2.getIndex() * CR_P.PARTY_REPLICATION, 
                                right_att_index * CR_P.PARTY_REPLICATION,
                                r_b, r_a, 
                                left_sum_index * CR_P.PARTY_REPLICATION);
                }
                output = new_r1;
                output.aSelectionUsed = true;
                free(r_a);
                free(r_b);
            } 
            else {
                leaderLogging("[Error] unsupported partial aggregation: " + agg->name + ".\n", true);
            }
        }
        else {

            auto sel = DerivedColumn<BShare, BShareTable>(output.getSize());
            auto sel_shares = sel.getShares();
            auto sel_shares_ptr = sel_shares.get();

            if (inputs[1].bSelectionUsed) {
                in_sel_right(&inputs[0].shareTable, &inputs[1].shareTable, c1.getIndex() * CR_P.PARTY_REPLICATION,
                            c2.getIndex() * CR_P.PARTY_REPLICATION,
                            inputs[1]["[SEL]"].getIndex() * CR_P.PARTY_REPLICATION,
                            sel_shares_ptr[0], inputs[0].getSize());
            } else {
                in(&inputs[0].shareTable, &inputs[1].shareTable, c1.getIndex() * CR_P.PARTY_REPLICATION,
                c2.getIndex() * CR_P.PARTY_REPLICATION,
                sel_shares_ptr[0], inputs[0].getSize());
            }

            // one communication round
            exchangeArrayWithParties<BShare>(sel_shares_ptr[0], sel_shares_ptr[1], output.getSize(), 1,
                                            PRIMITIVES_SHARE_TAG, 1);

            if (inputs[0].bSelectionUsed) {
            output["[SEL]"] = inputs[0]["[SEL]"] & sel;
            } else {
                output["[SEL]"] = sel;
            }
            
            output.bSelectionUsed = true;    
        }

        return true;
    }

public:
    std::vector<std::shared_ptr<Expression>> partialAgg;

    /**
     * The SemiJoinOperator updates the Boolean Selection Bit "[SEL]" for the left table
     *  based on an expression evaluation involving another table.
     * @param input_1 - The previous operator which this operator will process its output as left table input.
     * @param input_2 - The previous operator which this operator will process its output as right table input.
     * @param expression - An expression to specify the semi-join condition.
     * @param batchSize - Batch size this operator will use to process the input.
     */
    SemiJoinOperator(std::shared_ptr<Operator> input_1, std::shared_ptr<Operator> input_2,
                     const std::shared_ptr<Expression> expression,
                     int batchSize) : Operator(OperatorType::SemiJoin, batchSize) {
        operators.push_back(input_1);
        operators.push_back(input_2);
        expressions.push_back(expression);
    }


    /**
     * The SemiJoinOperator updates the Boolean Selection Bit "[SEL]" for the left table
     *  based on an expression evaluation involving another table. The operator also applies
     *  partial aggregations if any specified.
     * @param input_1 - The previous operator which this operator will process its output as left table input.
     * @param input_2 - The previous operator which this operator will process its output as right table input.
     * @param expression - An expression to specify the semi-join condition.
     * @param partialAgg - List of Expressions to specify additional aggregations.
     * @param batchSize - Batch size this operator will use to process the input (semi-join, aggregation).
     */
    SemiJoinOperator(std::shared_ptr<Operator> input_1, std::shared_ptr<Operator> input_2,
                     const std::shared_ptr<Expression> expression,
                     const std::vector<std::shared_ptr<Expression>> &partialAgg,
                     int batchSize) : Operator(OperatorType::SemiJoin, batchSize) {
        operators.push_back(input_1);
        operators.push_back(input_2);
        expressions.push_back(expression);
        this->partialAgg = partialAgg;
    }

};


#endif //SECRECY_SEMI_JOIN_OPERATOR_H