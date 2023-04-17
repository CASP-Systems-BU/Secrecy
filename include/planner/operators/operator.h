#ifndef SECRECY_OPERATOR_H
#define SECRECY_OPERATOR_H

#include <set>
#include <vector>
#include <iostream>

#include "expression.h"
#include "../relation.h"

#define B_SEL "[SEL]"
#define A_SEL "SEL"

enum OperatorType {
    NoOperator,
    Composition,
    Scan,
    Selection,
    Distinct,
    GroupBy, Join,
    SemiJoin,
    Sort,
    Conversion,
    Limit,
    Projection,
    TableFunction,
    Mask
};

static std::map<OperatorType, std::string> OperatorTypeName = {
        {OperatorType::NoOperator,  "NoOperator"},
        {OperatorType::Composition, "Composition"},
        {OperatorType::Scan,        "Scan"},
        {OperatorType::Selection,   "Selection"},
        {OperatorType::Distinct,    "Distinct"},
        {OperatorType::GroupBy,     "GroupBy"},
        {OperatorType::Join,        "Join"},
        {OperatorType::SemiJoin,    "SemiJoin"},
        {OperatorType::Sort,        "Sort"},
        {OperatorType::Conversion,  "Conversion"},
        {OperatorType::Limit,       "Limit"},
        {OperatorType::Projection,  "Projection"},
        {OperatorType::TableFunction,  "TableFunction"},
        {OperatorType::Mask,       "Mask"}
};


/**
 *
 */
class Operator : std::enable_shared_from_this<Operator>{
public:
    Relation output;
    friend class OptimizationRule;

    // Start: Should be private.
    OperatorType operatorType;
    int batchSize;

    std::vector<std::shared_ptr<Expression>> expressions;
    std::vector<std::shared_ptr<Operator>> operators;

    std::vector<Relation> inputs;
    // End: Should be private.

    /**
     * Starts executing operators in a DFS order.
     * @param level How many levels down from the root this operator is.
     * It is used organizing indentation in log printing.
     * @return Boolean to indicate execution is done successfully.
     */
    bool execute(int level = 0) {
        // execute children
        for (auto op : operators) {
            op->execute(level + 1);
            inputs.push_back(op->getOutput());
        }

        leaderLogging(getString(level));
        internalExecute();

        // // IMPORTANT NOTE: Uncomment this for debug logging
        // printOutputRelation(level);

        return true;
    }

    std::shared_ptr<Operator> getFirst(OperatorType opType, bool includeThis = true) {

        if (includeThis && this->operatorType == opType) {
            return shared_from_this();
        }

        for (auto op : this->operators) {
            if (op->operatorType == opType) {
                return op;
            } else {
                std::shared_ptr<Operator> child = op->getFirst(opType);
                if (child != nullptr) {
                    return child;
                }
            }
        }

        return nullptr;
    }

    void printOperator(int level = 0) {
        leaderLogging(getString(level));
        // execute children
        for (auto op : operators) {
            op->printOperator(level + 1);
        }
    }

    // Name of the relation this operator is attached to.
    std::string getName() const {
        switch (this->operatorType) {
            case OperatorType::Scan:
                return inputs[0].getName();
                break;
            case OperatorType::Join:
                return operators[0]->getName() + "*" + operators[1]->getName();
                break;
            default:
                return operators[0]->getName();
                break;
        }
    }

    // This operator affects one of the mentioned tables
    virtual bool hasTableDependency(const std::map<std::string, std::set<std::string>> &tables) {
        auto tableName = getName();

        for (auto it = tables.begin(); it!=tables.end(); ++it){
            if(tableName.find(it->first) != std::string::npos){
                return true;
            }
        }

        return false;
    }

    virtual bool hasTableDependency(const std::vector<std::shared_ptr<Expression>> &exprs) {
        auto tableName = getName();

        for (auto it = exprs.begin(); it!=exprs.end(); ++it){
            if(tableName.find((*it)->toString()) != std::string::npos){
                return true;
            }
        }

        return false;
    }

    // This operator affects one of the specific mentioned columns
    virtual bool hasColumnDependency(std::map<std::string, std::set<std::string>> &tables, bool notDirect = true) {
        switch (this->operatorType) {
            case OperatorType::Scan:
                return hasTableDependency(tables);
                break;
            default:
                for (auto expr : this->expressions) {
                    if (tables.count(expr->toString()) &&
                        tables[expr->toString()].count(expr->toString())) {
                        return true;
                    }
                }
        }

        if (notDirect) {
            // recursion
            for (auto op:this->operators) {
                if (op->hasColumnDependency(tables)) {
                    return true;
                }
            }
        }
        return false;
    }


    virtual bool hasColumnDependency(const std::vector<std::shared_ptr<Expression>> &exprs, bool notDirect = true) {
        switch (this->operatorType) {
            case OperatorType::Scan:
                return hasTableDependency(exprs);
                break;
            default:
                for (auto expr : this->expressions) {
                    for(auto expr2: exprs){
                        if (expr->toString() == expr2->toString()
                            && expr->toString() == expr2->toString()) {
                            return true;
                        }
                    }

                }
        }

        if (notDirect) {
            // recursion
            for (auto op:this->operators) {
                if (op->hasColumnDependency(exprs)) {
                    return true;
                }
            }
        }
        return false;
    }

    virtual Schema getOutputSchema() {
        std::cerr << "[ERROR]: getOutputSchema() not implemented." << std::endl;
        exit(-1);
    }

    Relation &getOutput() {
        return output;
    }


public:

    Operator() {
        this->operatorType = OperatorType::NoOperator;
        this->batchSize = -1;
    }

    Operator(const OperatorType &operatorType, const int &batchSize = 8) {
        this->operatorType = operatorType;
        this->batchSize = batchSize;
    }

    virtual ~Operator(){}


    virtual bool internalExecute() {
        return false;
    }

    virtual std::string getString(int indentNumber = 0, bool children = false) {
        std::string str = std::string(indentNumber * 2, ' ');
        str += OperatorTypeName[operatorType];

        // Print operator specific info
        str += "(" + getInternalString() + ")";

        // Print expressions
        str += "[";
        for (auto expression:expressions) {
            str += expression->toString() + ", ";
        }
        if (expressions.size()) {
            str = str.substr(0, str.size() - 2);
        }

        str += "]";

        str += "\n";
        return str;
    }

    virtual std::string getInternalString() {
        return "";
    }


    void printOutputRelation(const int& level) {
        // print table
        auto opened_ = output.openRelation();
        auto opened_ptr = opened_.get();
        if (get_rank() == 0) {

            // Header:
            auto cols = output.getCols();
            std::string header_str = std::string(level * 3, '\t');
            for (auto &col : cols) {
                header_str += ((col) + (output.isBColumnName(col) ? "(B)" : "(A)") + ("\t"));
            }
            header_str += "\n";
            leaderLogging(header_str);

            // Data:
            for (int i = 0; i < output.shareTable.numRows; ++i) {
                std::cout << std::string(level * 3, '\t');
                printf("%d:\t", i);
                for (int j = 0; j < output.shareTable.numCols; ++j) {
                    printf("%lld\t", opened_ptr[j][i]);
                }
                printf("\n");
            }
        }
    }

};


#endif //SECRECY_OPERATOR_H