#ifndef SECRECY_JOIN_OPERATOR_H
#define SECRECY_JOIN_OPERATOR_H

#include "operator.h"
#include "../../core/relational.h"
#include "../../planner/table/column_conversion.h"

/**
 * The JoinOperator merges two tables based on a merge expression.
 * The operators constructs new table based on the used batch sizes.
 * For each left table batch size, the whole right table will joined using
 * the right table batch size.
 * Note: If input tables are sorted and you want to keep this sorting,
 * use left table batch size of 1.
 */
class JoinOperator : public Operator {
    int batch_size_1;
    int batch_size_2;

    virtual bool internalExecute() {
        auto r1 = inputs[0];
        auto r2 = inputs[1];

        Column &c1 = r1[expressions[0]->expressions[0]->toString()];
        Column &c2 = r2[expressions[0]->expressions[1]->toString()];
        const int c_1 = c1.getIndex() * CR_P.PARTY_REPLICATION;
        const int c_2 = c2.getIndex() * CR_P.PARTY_REPLICATION;

        // Get New Size
        const int rows_1 = r1.getSize();
        const int rows_2 = r2.getSize();
        int rows = rows_1 * rows_2;

        // Get New Name
        std::string name = r1.getName() + "*" + r2.getName();

        // Merge column Names
        std::vector<std::string> cols = r1.getCols();
        int cols_1 = cols.size();
        for (auto col: r2.getCols()) {
            cols.push_back(col);
        }
        int cols_2 = cols.size() - cols_1;

        // One more selection for the large table
        cols.push_back("SEL");
        cols.push_back("[SEL]");

        // Create New relation
        output = Relation(name, cols, rows);
        output.allocateTable();

        int len_batch = batch_size_1 * batch_size_2;

        // Merge Tables
        int base_col = CR_P.PARTY_REPLICATION * cols_1;
        for (int i = 0; i < output.shareTable.numRows; ++i) {
            int cols_1_ = CR_P.PARTY_REPLICATION * cols_1;
            for (int j = 0; j < cols_1_; ++j) {
                int i_ = (i / batch_size_2) % batch_size_1 +
                         (i / (batch_size_1 * inputs[1].shareTable.numRows)) * batch_size_1;
                output.shareTable.content[i][j] = inputs[0].shareTable.content[i_][j];
            }
            int cols_2_ = CR_P.PARTY_REPLICATION * cols_2;
            for (int j = 0; j < cols_2_; ++j) {
                int i_ = (i % (batch_size_1 * inputs[1].shareTable.numRows));
                i_ = i_ % batch_size_2 + i_ / len_batch * batch_size_2;
                int j_ = base_col + j;
                output.shareTable.content[i][j_] = inputs[1].shareTable.content[i_][j];
            }
        }

        // For output "[SEL]"
        DerivedColumn<BShare, BShareTable> sel(rows);
        auto sel_column_array = sel.getShares();
        auto sel_column_array_ptr = sel_column_array.get();

        switch (expressions[0]->expressionType) {
            case Expression::Equal:

                for (int i = 0; i < rows_1; i += batch_size_1) {
                    for (int j = 0; j < rows_2; j += batch_size_2) {
                        int end_1 = i + batch_size_1 - 1 < rows_1 ? i + batch_size_1 - 1 : rows_1 - 1;
                        int end_2 = j + batch_size_2 - 1 < rows_2 ? j + batch_size_2 - 1 : rows_2 - 1;
                        int offset = i * rows_2 + j * batch_size_1;

                        join_eq_b_batch(&r1.shareTable, &r2.shareTable,
                                        i, end_1 + 1, j, end_2 + 1, c_1, c_2,
                                        &sel_column_array_ptr[1][offset], &sel_column_array_ptr[0][offset]);
                    }
                }

                break;
            case Expression::GreaterThanEqual:
                for (int i = 0; i < rows_1; i += batch_size_1) {
                    for (int j = 0; j < rows_2; j += batch_size_2) {
                        int end_1 = i + batch_size_1 - 1 < rows_1 ? i + batch_size_1 - 1 : rows_1 - 1;
                        int end_2 = j + batch_size_2 - 1 < rows_2 ? j + batch_size_2 - 1 : rows_2 - 1;
                        join_geq_b_batch(&r1.shareTable, &r2.shareTable,
                                         i, end_1, j, end_2, c_1, c_2,
                                         sel_column_array_ptr[1], sel_column_array_ptr[0]);
                    }
                }
                break;
            default:
                break;
        }

        // Exchange communication
        exchangeArrayWithParties<BShare>(sel_column_array_ptr[0], sel_column_array_ptr[1], rows, 1,
                                         PRIMITIVES_SHARE_TAG, 1);

        output.aSelectionUsed = inputs[0].aSelectionUsed | inputs[1].aSelectionUsed;
        if (inputs[0].aSelectionUsed & inputs[1].aSelectionUsed) {
            output["SEL"] = output[inputs[0].getName() + ".SEL"] * output[inputs[0].getName() + ".SEL"];
        } else if (inputs[0].aSelectionUsed) {
            output["SEL"] = output[inputs[0].getName() + ".SEL"];
        } else if (inputs[1].aSelectionUsed) {
            output["SEL"] = output[inputs[1].getName() + ".SEL"];
        }

        output.bSelectionUsed = true;
        if (inputs[0].bSelectionUsed & inputs[1].bSelectionUsed) {
            output["[SEL]"] = output[inputs[0].getName() + ".[SEL]"] & output[inputs[1].getName() + ".[SEL]"] & sel;
        } else if (inputs[0].bSelectionUsed) {
            output["[SEL]"] = output[inputs[0].getName() + ".[SEL]"] & sel;
        } else if (inputs[1].bSelectionUsed) {
            output["[SEL]"] = output[inputs[1].getName() + ".[SEL]"] & sel;
        } else {
            output["[SEL]"] = sel;
        }

        return true;
    }

public:

    /**
     * The JoinOperator merges two tables based on a merge expression.
     * The operators constructs new table based on the used batch sizes.
     * For each left table batch size, the whole right table will joined using
     * the right table batch size.
     * Note: If input tables are sorted and you want to keep this sorting,
     * use left table batch size of 1.
     * @param input_1 - The previous operator which this operator will process its output as left table input.
     * @param input_2 - The previous operator which this operator will process its output as right table input.
     * @param expression - An expression to specify the join condition.
     * @param batch_size_1 - Batch size this operator will use to process left input table. (set to 1 to keep it sorted.)
     * @param batch_size_2 - Batch size this operator will use to process right input table.
     */
    JoinOperator(std::shared_ptr<Operator> input_1, std::shared_ptr<Operator> input_2,
                 const std::shared_ptr<Expression> expression, int batch_size_1 = 8, int batch_size_2 = 8) :
            Operator(OperatorType::Join, batch_size_2) {
        operators.push_back(input_1);
        operators.push_back(input_2);
        expressions.push_back(expression);
        this->batch_size_1 = batch_size_1;
        this->batch_size_2 = batch_size_2;
    }


};

#endif //SECRECY_JOIN_OPERATOR_H