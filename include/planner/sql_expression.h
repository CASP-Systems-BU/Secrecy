#ifndef SECRECY_SQL_EXPRESSION_H
#define SECRECY_SQL_EXPRESSION_H

#include <string>

#include "costs.h"
#include "operators/expression.h"
#include "relation.h"

class Problem; // Forward declaration of Problem in the DP table

/***
 * A SQL expression (i.e., a logical operator)
 */
class SQLExpression : public std::enable_shared_from_this<SQLExpression> {
public:
    enum SQLExpressionType { Scan,
                            Join,
                            Filter,
                            GroupBy,
                            OrderBy,
                            Distinct,
                            TableFunction,
                            Unknown
    };
    /***
     * Expression id
     */
    int id;
    /***
     * Expression type
     */
    SQLExpressionType eType;
    /***
     * Constructor
     * @param id - The id of the expression (optional).
     */
    SQLExpression(int id=-1, SQLExpressionType eType=Unknown) : id(id), eType(eType) { };
    /***
     * Destructor
     */
    ~SQLExpression() {};

    /***
     * Returns a string of the expression (must be implemented by the sub-classes).
     */
    virtual std::string toString() const = 0;
};

class ScanExpr : public SQLExpression {
public:
    /***
     * The scanned relation
     */
    std::shared_ptr<RelationExpr> rel;

    /***
     * Constructor
     * @param id - The id of the expression (optional).
     */
    ScanExpr(int id=-1) : SQLExpression(id, Scan) { };
    /***
     * Destructor
     */
    ~ScanExpr() { };

    std::string toString() const;

    friend class PlanGen;
};

/***
 * A join expression
 */
class JoinExpr : public SQLExpression {
public:
    /***
     * Relational join type
     */
    enum JoinType {
        Inner,
        LeftOuter,
        RightOuter,
        Full,
        Semi
    };
    /***
     * Join types
     */
    const std::map<JoinType, std::string> JoinTypeName = {
            {JoinType::Inner,         "InnerJoin"},
            {JoinType::LeftOuter,     "LeftOuterJoin"},
            {JoinType::RightOuter,    "RightOuterJoin"},
            {JoinType::Full,          "FullOuterJoin"},
            {JoinType::Semi,          "SemiJoin"}
    };

    /***
     * The join type
     */
    JoinType jType;
    /***
     * The join predicate
     */
    std::shared_ptr<Expression> predicate;
    /***
     * Partial aggregations applied within a semi-join (if any)
     */
    std::vector<std::shared_ptr<Expression>> agg;
    /***
     * The left input
     */
    std::shared_ptr<RelationExpr> left;
    /***
     * The right input
     */
    std::shared_ptr<RelationExpr> right;
    /***
     * Left and right relation ids TODO: merge with left and right
     */
    int rId1 = -1, rId2 = -1;

    /***
     * Constructor
     * @param jType - The join type
     * @param predicate - The join predicate
     * @param left - The left input
     * @param right - The right input
     * @param id - The join id (optional)
     * @param rId1 - The left input id (optional)
     * @param rId2 - The right input id (optional)
     */
    JoinExpr(JoinType jType, std::shared_ptr<Expression> predicate, std::shared_ptr<RelationExpr> left,
             std::shared_ptr<RelationExpr> right, int id=-1, int rId1=-1, int rId2=-1)
                : SQLExpression(id, Join), jType(jType), predicate(predicate), left(left), right(right),
                  rId1(rId1), rId2(rId2) { };

    /***
     * Checks if the join predicate can be applied to create a larger plan.
     * @param p1 - The first problem in the DP table the join predicate applies to.
     * @param p2 - The second problem in the DP table the join predicate applies to.
     * @return True if the join applies, False otherwise.
     */
    bool applies(const std::shared_ptr<Problem> p1, const std::shared_ptr<Problem> p2) const;
    /***
     * Checks if the join predicate can be applied as a filter.
     * @param p - The problem in the DP table the join predicate applies to.
     * @return True if the join applies, False otherwise.
     */
    bool applies(const std::shared_ptr<Problem> p) const;
    /***
     * @return True if `this` is a semi-join, False otherwise.
     */
    bool semi() { return jType==JoinType::Semi; };
    /***
     * @return True if `this` is a self-join, False otherwise.
     */
    bool self() const;
    /***
     * Computes and returns the cost of the join.
     * @param cardinality - The cardinality of the cartesian product of the input relations (optional)
     * @return The cost of the join.
     */
    CostModel::Cost cost(CostModel::card_t cardinality=1) const {
        switch (this->predicate->expressionType) {
            case Expression::Equal:
            case Expression::In:
                return CostModel::eqCost(cardinality);
            case Expression::GreaterThan:
            case Expression::GreaterThanEqual:
            case Expression::LessThan:
            case Expression::LessThanEqual:
                return CostModel::ineqCost(cardinality);
            default:
                leaderLogging("[ERROR] Unsupported join predicate.\n", true);
        }
        // Shouldn't reach this point
        return CostModel::Cost();
    };
    /***
     * Returns the root of the tree the join expression belongs to.
     */
     std::shared_ptr<Expression> findRoot() const;
    /***
     * Prints join info to stdout.
     */
    void print() const;
    /***
     * @return Join expression as a string.
     */
    std::string toString() const override;
};

/***
 * A filter expression
 */
class FilterExpr : public SQLExpression {
public:
    /***
     * The type
     */
    Expression::ExpressionType type;
    /***
     * Attributes the filter is applied to.
     */
    std::vector<std::shared_ptr<Expression>> attributes;
    /***
     * The actual predicate
     */
    std::shared_ptr<Expression> selPredicate;
    /***
     * The id of the relation the filter is applied to.
     */
    unsigned rId;

    /***
     * Constructor
     * @param type - The filter type
     * @param attributes - The attributes the filter is applied to.
     * @param predicate - The filter predicate.
     */
    FilterExpr(Expression::ExpressionType type,
               std::vector<std::shared_ptr<Expression>>& attributes,
               std::shared_ptr<Expression> predicate)
            : SQLExpression(-1, Filter), type(type), attributes(attributes), selPredicate(predicate) { }
    /***
     * Shallow copy constructor
     * @param other - The filter whose contents will be used in the new filter.
     */
    FilterExpr(const FilterExpr &other)
            : SQLExpression(-1, Filter), type(other.type), attributes(other.attributes),
              selPredicate(other.selPredicate) { };

    /***
     * Destructor
     */
    ~FilterExpr() { };

    /***
     * Filter assignment
     * @param other - The filter to refer to.
     * @return
     */
    FilterExpr& operator=(const FilterExpr &other) { return *this; };

    /***
     * Checks if the filter can be applied to create a larger plan.
     * @param prob - The problem in the DP table the filter can be applied to.
     * @return - True if the filter can be applied, False otherwise.
     */
    bool applies(const std::shared_ptr<Problem> prob) const;
    /***
     * Computes and returns the cost of the filter.
     * @param cardinality - The cardinality of the input relation (optional)/
     * @return The cost of the filter.
     */
    CostModel::Cost cost(CostModel::card_t cardinality=1) const {
        switch (this->type) {
            case Expression::Equal:
                return CostModel::eqCost(cardinality);
            case Expression::GreaterThan:
            case Expression::GreaterThanEqual:
            case Expression::LessThan:
            case Expression::LessThanEqual:
                return CostModel::ineqCost(cardinality);
            default:
                leaderLogging("[ERROR] Unsupported filter.\n", true);
        }
        // Shouldn't reach this point
        return CostModel::Cost();
    }
    /***
     * Prints filter info to stdout.
     */
    void print() const;
    /***
     * @return The filter expression as a string.
     */
    std::string toString() const override;
};

/***
 * A group-by expression
 */
class GroupByExpr : public SQLExpression {
public:
    /***
     * Group-by keys
     */
    std::vector<std::shared_ptr<Expression>> keys;
    /***
     * Aggregation expressions
     */
    std::vector<std::shared_ptr<Expression>> aggs;
    /***
     * Having expression
     */
    std::shared_ptr<Expression> having;
    /***
     * True if `this` is an odd-even aggregation, False otherwise.
     */
    bool oddEven = false;

    /***
     * Constructor
     */
    GroupByExpr(int id=-1) : SQLExpression(id, GroupBy) { };
    /***
     * Destructor
     */
    ~GroupByExpr() {};

    /***
     * @return The group-by expression as a string.
     */
    std::string toString() const override;
};

/***
 * An distinct expression
 */
class DistinctExpr : public SQLExpression {
public:
    /***
     * Distinct keys
     */
    std::vector<std::shared_ptr<Expression>> keys;
    /***
     * True if `this` is an AdjacentEq expression, False otherwise
     */
    bool adjacentEq = false;

    /***
     * Constructor
     * @param id - The expression id.
     */
    DistinctExpr(int id=-1) : SQLExpression(id, Distinct) { };
    /***
     * Destructor
     */
    ~DistinctExpr() { };

    /***
     * @return The distinct expression as a string.
     */
    std::string toString() const override;
};

/***
 * An order-by expression
 */
class OrderByExpr : public SQLExpression {
public:
    /***
     * Sorting keys
     */
    std::vector<std::shared_ptr<Expression>> keys;
    /***
     * Sorting direction per key (True for ASC, False for DESC)
     */
    std::vector<bool> asc;
    /***
     * The LIMIT
     */
    int limit;

    /***
     * Constructor
     * @param id - The expression id
     */
     OrderByExpr(int id=-1) : SQLExpression(id, OrderBy) { };
    /***
     * Destructor
     */
    ~OrderByExpr() { };

    /***
     * @return The order-by expression as a string
     */
    std::string toString() const override;
};

/***
 * A table function expression
 */
class TableFunctionExpr : public SQLExpression {
public:
    /***
     * The function expression
     */
    std::shared_ptr<Expression> exp;
    /***
     * True if `this` is a masking function, False otherwise
     */
    bool mask = false;
    /***
     * True if the expression is applied to adjacent records, False otherwise
     */
    bool adjacent = false;

    /***
     * Constructor
     * @param exp - The function expression
     */
    TableFunctionExpr(std::shared_ptr<Expression> exp, int id=-1) : SQLExpression(id, TableFunction), exp(exp) { };

    /***
     * Shallow copy constructor
     * @param other - The function expression to copy from.
     */
    TableFunctionExpr(const std::shared_ptr<TableFunctionExpr> other) : SQLExpression(-1, TableFunction) {
        id = other->id;
        mask = other->mask;
        adjacent = other->adjacent;
        auto p = std::dynamic_pointer_cast<FunctionExpr>(other->exp);
        assert(p);
        exp = std::make_shared<FunctionExpr>(FunctionExpr(p));
    };

    /***
     * Destructor
     */
    ~TableFunctionExpr() { };

    /***
     * Computes and returns the cost of the function.
     * @param cardinality - The cardinality of the input (optional).
     * @return The cost of the table function.
     */
    CostModel::Cost cost(CostModel::card_t cardinality=1) { return CostModel::Cost(); };
    /***
     * @return The function expression as a string.
     */
    std::string toString() const override;
};

#endif //SECRECY_SQL_EXPRESSION_H
