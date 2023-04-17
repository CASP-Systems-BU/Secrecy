#ifndef SECRECY_MASK_OPERATOR_H
#define SECRECY_MASK_OPERATOR_H

#include "operator.h"
#include "../../core/relational.h"

/**
 * The MaskOperator uses the Boolean Selection Bit "[SEL]"
 * to decide to either keep a row's values or replace them by ones (masking).
 * Note: assumes existence of valid column "[SEL]".
 * Note: operator modifies input table in place.
 */
class MaskOperator : public Operator {

    virtual bool internalExecute() {
        output = inputs[0];

        DerivedColumn<BShare, BShareTable> sel(output.getSize());
        sel = output["[SEL]"];
        auto sel_shares = sel.getShares();
        auto sel_shares_ptr = sel_shares.get();

        output.shareTable.numCols *= CR_P.PARTY_REPLICATION;
        mask(&output.shareTable, sel_shares_ptr[0], batchSize);
        output.shareTable.numCols /= CR_P.PARTY_REPLICATION;
        return true;
    }

    virtual std::string getInternalString() {
        return inputs[0].getName();
    }

public:

    /**
     * The MaskOperator uses the Boolean Selection Bit "[SEL]"
     * to decide to either keep a row's values or replace them by ones (masking).
     * Note: assumes existence of valid column "[SEL]".
     * Note: operator modifies input table in place.
     * @param input The previous operator which this operator will process its output.
     * @param batchSize Batch size this operator will use to process input data.
     */
    MaskOperator(std::shared_ptr<Operator> input, int batchSize) :
            Operator(OperatorType::Mask, batchSize) {
        this->operators.push_back(input);
    }

};


#endif //SECRECY_MASK_OPERATOR_H