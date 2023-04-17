#ifndef SECRECY_TABLE_FUNCTION_H
#define SECRECY_TABLE_FUNCTION_H

#include "operator.h"
#include "../../core/relational.h"

enum FunctionType {
    CNT, MX, MN, SM, AV, MN_MX, MUL, CUSTOM
};

static std::map<FunctionType, std::string> FunctionTypeName = {
        {FunctionType::CNT,        "COUNT"},
        {FunctionType::MX,          "MAX"},
        {FunctionType::MN,          "MIN"},
        {FunctionType::SM,          "SUM"},
        {FunctionType::AV,          "AVG"},
        {FunctionType::MUL,         "MUL"},
        {FunctionType::CUSTOM,       "CUSTOM"}
};

static std::map<std::string, FunctionType> FunctionNameType = {
        {"COUNT",       FunctionType::CNT},
        {"MAX",         FunctionType::MX},
        {"MIN",         FunctionType::MN},
        {"SUM",         FunctionType::SM},
        {"AVG",         FunctionType::AV},
        {"MUL",         FunctionType::MUL},
        {"CUSTOM",    FunctionType::CUSTOM}
};


class TableFunctionOperator : public Operator {

    virtual bool internalExecute() {
        // Pass output relation
        output = inputs[0];

        for(auto& columnName : inputs[0].getCols(true)){
            std::cout << columnName << "\t";
        }
        std::cout << std::endl;


        switch (functionType) {
            case FunctionType::CNT:
            {
                if(!output.aSelectionUsed){
                    output["SEL"] = convertToAColumn<AShare, AShareTable>(output["[SEL]"]);
                    output.aSelectionUsed = true;
                }
                auto sel_a = reinterpret_cast<const AColumn<AShare, AShareTable> *>(&output["SEL"]);

                output = Relation(output.getName(), {"SEL"}, 1);
                output["SEL"] = sel_a->sum();

                break;
            }
            case FunctionType::SM:
            {
                auto attr_a = reinterpret_cast<const AColumn<AShare, AShareTable> *>(&output[attributes[0]]);
                output = Relation(output.getName(), {"SUM"}, 1);
                output["SUM"] = attr_a->sum();

                break;
            }
            case FunctionType::MUL:
            {
                //TODO: assuming that the first expression is the result and then the two operands.
                output[attributes[0]] = output[attributes[1]] * output[attributes[2]];
                break;
            }
            default:{
                leaderLogging("[ERROR] Unsupported table function.\n", true);
            }

        }
        return true;
    }

    virtual std::string getInternalString() {
        return "TableFunction:" + std::to_string(cmpBit);
    }

public:
    FunctionType functionType;
    std::shared_ptr<Expression> functionExpr;
    std::vector<std::shared_ptr<Expression>> attributes;
    bool adjacent;
    bool cmpBit;

    TableFunctionOperator(std::shared_ptr<Operator> input, const FunctionType &functionType,
                          const std::shared_ptr<Expression> function_exp,
                          bool cmpBit=false, bool adjacent=false, int batchSize=8)
                          : Operator(OperatorType::TableFunction, batchSize)
    {
        this->operators.push_back(input);
        this->functionType = functionType;
        this->functionExpr = function_exp;
        this->cmpBit = cmpBit;
        this->adjacent = adjacent;
        this->batchSize = batchSize;
    };

    virtual std::string getString(int indentNumber = 0, bool children = false) {
        return "";
    }
};

#endif //SECRECY_TABLE_FUNCTION_H
