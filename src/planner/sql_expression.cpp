#include "../../include/planner/sql_expression.h"
#include "../../include/planner/dpTable.h"

// Checks if join predicate can be applied
bool JoinExpr::applies(const std::shared_ptr<Problem> p1,
             const std::shared_ptr<Problem> p2) const {
    if ((p1->relations.contains(this->rId1) && p2->relations.contains(this->rId2)) ||
        (p2->relations.contains(this->rId1) && p1->relations.contains(this->rId2)))
        return true;
    return false;
}

// Checks if join predicate can be applied as a filter
bool JoinExpr::applies(const std::shared_ptr<Problem> p) const {
    if ((p->relations.contains(this->rId1) && p->relations.contains(this->rId2)))
        return true;
    return false;
}

// Checks if filter predicate can be applied
bool FilterExpr::applies(const std::shared_ptr<Problem> prob) const {
    if (prob->relations.contains(this->rId))
        return true;
    return false;
}

std::string ScanExpr::toString() const {
    return "Scan (ID: " + std::to_string(id) + ", Relation: " + rel->toString() + ")";
}

// Returns the root of the tree the join expression belongs to
std::shared_ptr<Expression> JoinExpr::findRoot() const {
    assert(predicate != nullptr);
    return predicate->findRoot();
}

// Returns join expression as a string
std::string JoinExpr::toString() const {
    std::string jt(JoinTypeName.find(this->jType)->second);
    std::string s = jt + " (ID: " + std::to_string(this->id) + ", Predicate: " + this->predicate->toString();
    if (!this->agg.empty()) { // There exists a partial aggregation
        s += ", Partial Aggregation: " + this->agg[0]->toString();
    }
    s += ")";
    return s;
}

// Prints join info to stdout
void JoinExpr::print() const {
    std::cout << "Join type: (";
    std::string s(JoinTypeName.find(this->jType)->second);
    std::cout << s << ")   Predicate: " << this->predicate->toString() << std::endl;
}

// Returns filter expression as a string
std::string FilterExpr::toString() const {
    return "Filter (ID: " + std::to_string(this->id) + ", Predicate: " + selPredicate->toString() + ")";
}

// Prints filter info to stdout
void FilterExpr::print() const {
    std::cout << "Filter attributes: ";
    std::string s = "(";
    for (auto a : this->attributes)
        s += a->toString() + ", ";
    s = s.substr(0, s.size() - 2);
    std::cout << s << ")   " << "Predicate: " << selPredicate->toString() << std::endl;
}

// Returns order-by expression as a string
std::string DistinctExpr::toString() const {
    std::string t = adjacentEq ? "AdjacentEq" : "Distinct";
    std::string s =  t + " (ID: " + std::to_string(this->id) + ", Keys: ";
    for (int i=0; i<keys.size(); i++) {
        s += keys[i]->toString() + ", ";
    }
    s = s.substr(0, s.length() - 2) + ")";
    return s;
}

// Returns order-by expression as a string
std::string OrderByExpr::toString() const {
    assert(keys.size() == asc.size());
    std::string s = "OrderBy (ID: " + std::to_string(this->id) + ", Keys: ";
    for (int i=0; i<keys.size(); i++) {
        std::string d = asc[i] ? "ASC" : "DESC";
        s += keys[i]->toString() + " (" + d + "), ";
    }
    s = s.substr(0, s.length() - 2) + ")";
    return s;
}

// Returns group-by expression as a string
std::string GroupByExpr::toString() const {
    std::string t = oddEven ? "OddEvenAgg" : "GroupBy";
    std::string s = t + " (ID: " + std::to_string(this->id) + ", Keys: ";
    for (auto k: keys)
        s += k->toString() + ", ";
    if (!aggs.empty()) {
        s += "Aggregations: ";
        for (auto a: aggs)
            s += a->toString() + ", ";
    }
    if (this->having) {
        s += "Having: " + this->having->toString() + ", ";
    }
    s = s.substr(0, s.length() - 2) + ")";
    return s;
}


// Returns function expression as a string
std::string TableFunctionExpr::toString() const {
    std::string s = mask ? "Mask" : "TableFunction";
    return s + " (ID: " + std::to_string(this->id) + ", " + this->exp->toString() + ")";
}