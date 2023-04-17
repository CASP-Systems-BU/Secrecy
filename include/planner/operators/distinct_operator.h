#ifndef SECRECY_DISTINCT_OPERATOR_H
#define SECRECY_DISTINCT_OPERATOR_H

#include "operator.h"
#include "../../core/wrappers/distinct.h"
#include "../../core/relational.h"

/**
 * The DistinctOperator Modifies the padded Boolean Selection Bit "[SEL]"
 * to reflect selection of the first occurrence of a unique tuple
 * in the mentioned target columns.
 * - In case the selection flag is raised in the input,
 * the selection flag is taken into consideration
 * and then updated to reflect the distinct result.
 * - Operators expected input should be sorted or sorting flag set raised.
 * Note: Operator does not track the Arithmetic Boolean Selection Bit "SEL"
 * Note: If sorting is specified, the Arithmetic columns will be corrupted.
 */
class DistinctOperator : public Operator {
    virtual bool internalExecute() {
        // pass info to output
        output = inputs[0];

        unsigned num_keys = expressions.size();

        // if the input operator has generated a selected bit
        // we need to take it into account when sorting
        
        int sel_pos = -1;
        auto selected = inputs[0].bSelectionUsed;
        if (selected) {
            num_keys++;
            sel_pos = output["[SEL]"].getIndex();
        }

        int offset = 0; // strating offset of non-selected-bit keys
        unsigned key_pos[num_keys];
        if (selected) {
            key_pos[0] = sel_pos * CR_P.PARTY_REPLICATION;
            offset = 1;
        }
        // get the position of the other keys
        for (int i=offset; i<num_keys; i++) {
            auto rel_name = expressions[i-offset]->relation;
            auto key_name = expressions[i-offset]->name;
            key_pos[i] = (output[rel_name + "." + key_name].getIndex()) * CR_P.PARTY_REPLICATION;
        }
        
        if (includeSorting) {
            bool asc[num_keys];
            for (int i = 0; i < num_keys; ++i) {
                asc[i] = 1;
            }
            if (batchSize == output.getSize()) {
                batchSize /= 2;
            }
            output.shareTable.numCols *= CR_P.PARTY_REPLICATION;
            bitonic_sort_batch(&output.shareTable, key_pos, num_keys, asc, batchSize);
            output.shareTable.numCols /= CR_P.PARTY_REPLICATION;
        }

        // This function has one communication round in the end
        auto sel = distinct_array<BShare, BShareTable>(&output.shareTable,
                                                       key_pos, num_keys,
                                                       batchSize - 1);
        auto sel_c = DerivedColumn<BShare, BShareTable>(sel, output.getSize(),
                                                        ColumnType::ARRAY_REFERENCE);

        if (inputs[0].bSelectionUsed) {
            output["[SEL]"] = output["[SEL]"] & sel_c;
        } else {
            output["[SEL]"] = sel_c;
        }
        output.bSelectionUsed = true;

        return true;
    }

    virtual std::string getInternalString() {
        return "Sort:" + std::to_string(includeSorting);
    }

public:
    bool includeSorting;

    /**
     * The DistinctOperator Modifies the padded Boolean Selection Bit "[SEL]"
     * to reflect selection of the first occurrence of a unique tuple
     * in the mentioned target columns.
     * - In case the selection flag is raised in the input,
     * the selection flag is taken into consideration
     * and then updated to reflect the distinct result.
     * - Operators expected input should be sorted or sorting flag set raised.
     * Note: Operator does not track the Arithmetic Boolean Selection Bit "SEL"
     * Note: If sorting is specified, the Arithmetic columns will be corrupted.
     * @param input_1 The previous operator which this operator will process its output.
     * @param expressions Expressions to specify the target columns to run distinct on.
     * @param includeSorting If set to true, input data will be sorted using expressions in ascending order.
     * @param batchSize Batch size this operator will use to process input data (distinct + sorting).
     */
    DistinctOperator(std::shared_ptr<Operator> input_1, std::vector<std::shared_ptr<Expression>> expressions, bool includeSorting = true,
                     int batchSize = 8) : Operator(OperatorType::Distinct, batchSize) {
        this->operators.push_back(input_1);
        this->expressions = expressions;
        this->includeSorting = includeSorting;
    }


};


#endif //SECRECY_DISTINCT_OPERATOR_H