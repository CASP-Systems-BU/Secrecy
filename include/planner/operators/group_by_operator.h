#ifndef SECRECY_GROUP_BY_OPERATOR_H
#define SECRECY_GROUP_BY_OPERATOR_H

#include "operator.h"
#include "../../core/relational.h"

enum GroupByType {
    COUNT, MAX, MIN, SUM, AVG, MIN_MAX, MULTIPLE
};

static std::map<GroupByType, std::string> GroupByTypeName = {
        {GroupByType::COUNT,        "COUNT"},
        {GroupByType::MAX,          "MAX"},
        {GroupByType::MIN,          "Min"},
        {GroupByType::MIN_MAX,      "MIN_MAX"},  // TODO (john): Treat this as MULTIPLE
        {GroupByType::SUM,          "SUM"},
        {GroupByType::AVG,          "AVG"},
        {GroupByType::MULTIPLE,     "MULTIPLE"}
};

static std::map<std::string, GroupByType> GroupByNameType = {
        {"COUNT",       GroupByType::COUNT},
        {"MAX",         GroupByType::MAX},         
        {"MIN",         GroupByType::MIN},
        {"MIN_MAX",     GroupByType::MIN_MAX},  // TODO (john): Treat this as MULTIPLE
        {"SUM",         GroupByType::SUM},
        {"AVG",         GroupByType::AVG},
        {"MULTIPLE",    GroupByType::MULTIPLE}
};

/**
 * The Group By Operator applies one the aggregation function to
 * the columns of the input relation. This operator can also be used
 * for the first and second phase of the join-group by aggregation.
 */
class GroupByOperator : public Operator {

    virtual bool internalExecute() {
        // pass output relation
        output = inputs[0];

        //auto indices_pair = output.getIndices(expressions);
        unsigned num_keys = expressions.size();
        unsigned key_pos[num_keys];

        // get key positions for sorting and grouping (must include replication)
        for (int i = 0; i < num_keys; i++) {
            auto rel_name = expressions[i]->relation;
            auto key_name = expressions[i]->name;
            key_pos[i] = (output[rel_name + "." + key_name].getIndex()) * CR_P.PARTY_REPLICATION;
        }

        if (includeSorting) {
            // sorting indices
            bool asc[num_keys];
            for (int i = 0; i < num_keys; ++i) {
                asc[i] = 1;
            }

            output.shareTable.numCols *= CR_P.PARTY_REPLICATION;
            bitonic_sort_batch(&output.shareTable, key_pos, num_keys, asc, batchSize);
            output.shareTable.numCols /= CR_P.PARTY_REPLICATION;
        }

        int r_size = output.getSize() * log2(output.getSize()) + 1;
        AShare *r_a = new AShare[r_size];
        BShare *r_b = new BShare[r_size];
        get_next_rb_array(r_b, r_size);
        get_next_array(r_a, r_size);

        switch (groupByType) {
            case GroupByType::COUNT: {
                if (output.aSelectionUsed && output.bSelectionUsed) {
                    DerivedColumn<AShare, AShareTable> c_sel_a(output.getSize());
                    c_sel_a = output["SEL"];
                    auto c_sel_a_shares = c_sel_a.getShares();
                    AShares *c_sel_a_shares_ptr = c_sel_a_shares.get();

                    DerivedColumn<BShare, BShareTable> c_sel_b(output.getSize());
                    c_sel_b = output["[SEL]"];
                    auto c_sel_b_shares = c_sel_b.getShares();
                    auto c_sel_b_shares_ptr = c_sel_b_shares.get();

                    output.shareTable.numCols *= CR_P.PARTY_REPLICATION;
                    group_by_count_sel_odd_even(&output.shareTable, key_pos, num_keys,
                                                batchSize, c_sel_b_shares_ptr[0], c_sel_a_shares_ptr[0],
                                                r_b, r_a);
                    output.shareTable.numCols /= CR_P.PARTY_REPLICATION;

                    // retrieve remote result
                    exchangeArrayWithParties<AShare>(c_sel_a_shares_ptr[0], c_sel_a_shares_ptr[1], output.getSize(),
                                                     1, PRIMITIVES_SHARE_TAG, 1);
                    exchangeArrayWithParties<BShare>(c_sel_b_shares_ptr[0], c_sel_b_shares_ptr[1], output.getSize(),
                                                     1, PRIMITIVES_SHARE_TAG, 1);

                    output["SEL"] = c_sel_a;
                    output["[SEL]"] = c_sel_b;
                } else {
                    leaderLogging(
                            "[ERROR] Unsupported use of group-by count: (aSelectionUsed && bSelectionUsed) must be true\n: ",
                            true);
                }
                break;
            }
            case GroupByType::SUM: {
                if (secondPhase) {
                    unsigned int sum_res_index = output["SEL"].getIndex() * CR_P.PARTY_REPLICATION;
                    // this is the second phase of aggregation decomposition
                    output.shareTable.numCols *= CR_P.PARTY_REPLICATION;
                    group_by_sum_odd_even(&output.shareTable, batchSize, r_b, r_a,
                                          sum_res_index, 0, key_pos, num_keys);
                    output.shareTable.numCols /= CR_P.PARTY_REPLICATION;
                } else {
                    leaderLogging("Unsupported use of group-by sum: secondPhase is false\n: ", true);
                }
                break;
            }
            default:
                break;
        }
        free(r_a);
        free(r_b);

        return true;
    }

    virtual std::string getInternalString() {
        return "Sort:" + std::to_string(includeSorting);
    }

public:
    GroupByType groupByType;
    std::vector<std::shared_ptr<Expression>> result_expressions;
    std::vector<GroupByType> aggFunctions;
    std::vector<std::shared_ptr<Expression>> aggAttributes;
    bool includeSorting;
    bool includeJoin;
    bool secondPhase;

    Expression joinExpression;

    /**
     * The Group By Operator applies one the aggregation function to
     * the columns of the input relation.
     * @param input_1 - The previous operator which this operator will process its output.
     * @param groupByType - The aggregation type (i.e. sum, max, min, count ...)
     * @param keyAttributes - A vector of expressions to specify the columns to apply the aggregation on.
     * @param resultAttributes - A vector of expressions to specify the columns to store the aggregation result in.
     * @param includeSorting - Whether to sort using the key attributes before applying the odd-even aggregation.
     * @param batchSize - Batch size this operator will use to process its input.
     */
    GroupByOperator(std::shared_ptr<Operator> input_1, const GroupByType &groupByType,
                    const std::vector<std::shared_ptr<Expression>> &keyAttributes,
                    const std::vector<std::shared_ptr<Expression>> &resultAttributes,
                    bool includeSorting = true, int batchSize = 8) :
            Operator(OperatorType::GroupBy, batchSize),
            joinExpression(Expression::Unknown) {
        this->operators.push_back(input_1);
        this->expressions = keyAttributes;
        this->result_expressions = resultAttributes;
        this->groupByType = groupByType;
        this->includeSorting = includeSorting;
        this->includeJoin = false;
    }

    /**
     * The Group By Operator applies one the aggregation function to
     * the columns of the input relation. This is the constructor for the second phase
     * of the join-group by aggregation.
     * @param input_1 - The previous operator which this operator will process its output.
     * @param groupByType - The aggregation type (i.e. sum, max, min, count ...)
     * @param keyAttributes - A vector of expressions to specify the columns to apply the aggregation on.
     * @param includeSorting - Whether to sort using the key attributes before applying the odd-even aggregation.
     * @param batchSize - Batch size this operator will use to process its input.
     * @param isSecondPhase Whether to apply the second phase of the join-group by aggregation.
     */
    GroupByOperator(std::shared_ptr<Operator> input_1, const GroupByType &groupByType,
                    const std::vector<std::shared_ptr<Expression>> &keyAttributes,
                    bool includeSorting = true, int batchSize = 8, bool isSecondPhase = false) :
            Operator(OperatorType::GroupBy, batchSize),
            joinExpression(Expression::Unknown) {
        this->operators.push_back(input_1);
        this->expressions = keyAttributes;
        this->groupByType = groupByType;
        this->includeSorting = includeSorting;
        this->includeJoin = false;
        this->secondPhase = isSecondPhase;
    }

     /**
      * The Group By Operator applies one the aggregation function to
      * the columns of the input relation. This is the constructor for the first phase
      * of the join-group by aggregation.
      * @param input_1 - The previous operator which this operator will process its output as left table input.
      * @param input_2 - The previous operator which this operator will process its output as right table input.
      * @param joinExpression - An expression to specify the join condition.
      * @param groupByType - The aggregation type (i.e. sum, max, min, count ...)
      * @param keyAttributes - A vector of expressions to specify the columns to apply the aggregation on.
      * @param resultAttributes - A vector of expressions to specify the columns to store the aggregation result in.
      * @param includeSorting - Whether to sort using the key attributes before applying the odd-even aggregation.
      * @param batchSize - Batch size this operator will use to process its input.
      */
    GroupByOperator(std::shared_ptr<Operator> input_1, std::shared_ptr<Operator> input_2,
                    const Expression &joinExpression,
                    const GroupByType &groupByType,
                    const std::vector<std::shared_ptr<Expression>> &keyAttributes,        // Has to be boolean shares
                    const std::vector<std::shared_ptr<Expression>> &resultAttributes,   // Can be either boolean shares or arithmetic shares
                    bool includeSorting = true, int batchSize = 8) :
            Operator(OperatorType::GroupBy, batchSize),
            joinExpression(joinExpression) {
        this->operators.push_back(input_1);
        this->expressions = keyAttributes;
        this->result_expressions = resultAttributes;
        this->groupByType = groupByType;
        this->includeSorting = includeSorting;

        this->includeJoin = true;
        this->operators.push_back(input_2);
    }

    virtual std::string getString(int indentNumber = 0, bool children = false) {
        std::string str = std::string(indentNumber * 2, ' ');
        str += OperatorTypeName[operatorType];

        // Print operator specific info
        str += "(" + getInternalString() + ")";

        // Print group-by keys
        str += "[(";
        for (auto expression:expressions) {
            str += expression->toString() + ", ";
        }
        if (expressions.size()) {
            str = str.substr(0, str.size() - 2);
        }
        str += "), ";
        // Print aggregations
        for (int i=0; i<aggFunctions.size(); i++) {
            str += GroupByTypeName[aggFunctions[i]] + "(" + aggAttributes[i]->toString() + "), ";
        }
        if (aggFunctions.size()) {
            str = str.substr(0, str.size() - 2);
        }
        str += "]\n";
        return str;
    }
};


#endif //SECRECY_GROUP_BY_OPERATOR_H