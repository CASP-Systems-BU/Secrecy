#ifndef SECRECY_PARSER_UTILS_H
#define SECRECY_PARSER_UTILS_H

#include <climits>

#include "../database.h"
#include "../operators/operators.h"
#include "../../external-lib/sql-parser/src/SQLParser.h"
#include "../../external-lib/sql-parser/src/util/sqlhelper.h"

class Expression;
class QueryGraph;

/***
 * Relational join type.
 */
enum JoinType {
    Inner,
    LeftOuter,
    RightOuter,
    Full,
    LeftSemi,
    RightSemi
};

/***
 * Join description.
 */
struct JoinDesc {
    /***
     * Join type
     */
    JoinType jType;
    /***
     * Join predicate
     */
    std::shared_ptr<Expression> predicate;
    /***
     * Left input
     */
    std::shared_ptr<RelationExpr> left;
    /***
     * Right input
     */
    std::shared_ptr<RelationExpr> right;
};

/***
 * Group-by description.
 */
struct GroupByDesc {
    /***
     * Group-by keys.
     */
    std::vector<std::shared_ptr<Expression>> keys;
    /***
     * Aggregation expressions.
     */
    std::vector<std::shared_ptr<Expression>> aggs;
    /***
     * Having expression.
     */
    std::shared_ptr<Expression> having;
};

/***
 * Order-by description.
 */
struct OrderByDesc {
    /***
     * List of sorting keys
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
};

/***
 * Query description.
 */
struct QueryDesc {
    /***
     * Query ID
     */
    int qId;
    /***
     * Parent query ID
     */
    int pId;
    /***
     * The list of attributes and expressions in the SELECT clause
     */
    std::vector<std::shared_ptr<Expression>> projections;
    /***
     * SELECT DISTINCT flag
     */
    bool distinct;
    /***
     * The list of relations involved in the query
     */
    std::vector<std::shared_ptr<RelationExpr>> relations;
    /***
     * The expression tree in the WHERE clause
     */
    std::shared_ptr<Expression> whereClause;
    /***
     * The list of filters in the query
     */
    std::vector<std::shared_ptr<Expression>> filters;
    /***
     * The list of joins in the query
     */
    std::vector<std::shared_ptr<JoinDesc>> joins;
    /***
     * Group-by clause
     */
    std::shared_ptr<GroupByDesc> groupBy;
    /***
     * Order-by clause
     */
    std::shared_ptr<OrderByDesc> orderBy;
};

/***
 * A parsed SQL query.
 */
struct SQLQuery {
    /**
     * The query id
     */
    int id;
    /**
     * The query string in SQL syntax.
     */
    std::string qString;
    /**
     * A pointer to the Abstract Syntax Tree (AST).
     */
    const hsql::SelectStatement *tree;
    /**
     * Pointers to ASTs of nested queries.
     */
    std::vector<const hsql::SelectStatement*> subTrees;
    /**
     * Returns a pointer to the SELECT clause of the query (if any).
     * @return A pointer to the AST of the SELECT clause
     */
    const std::vector<hsql::Expr*>* selectClause() const { return this->tree->selectList; };
    /**
     * @return True if there is a SELECT DISTINCT clause, False otherwise.
     */
    bool hasDistinct() const { return this->tree->selectDistinct; };
    /**
     * Returns a pointer to the FROM clause of the query.
     * @return A pointer to the AST of the FROM clause
     */
    const hsql::TableRef* fromClause() const { return this->tree->fromTable; };
    /**
     * Returns a pointer to the WHERE clause of the query.
     * @return A pointer to the AST of the WHERE clause
     */
    const hsql::Expr* whereClause() const { return this->tree->whereClause; };
    /**
     * Returns a pointer to the GROUP-BY clause of the query.
     * @return A pointer to the AST of the GROUP-BY clause
     */
    const hsql::GroupByDescription* groupByClause() const { return this->tree->groupBy; };
    /**
     * Returns a pointer to the ORDER-BY clause of the query.
     * @return A pointer to the AST of the ORDER-BY clause
     */
    const std::vector<hsql::OrderDescription*>* orderByClause() const { return this->tree->order; };
};

/***
 * Collects a parsed SQL expression and constructs the corresponding Secrecy expression.
 * @param expr - A pointer to the AST of the parsed expression.
 * @param query - The SQL query the expression belongs to.
 * @return The respective Secrecy expression.
 */
static std::shared_ptr<Expression> processExpression(const hsql::Expr *expr, SQLQuery& query) {
    assert(expr != nullptr);
    ExpressionType eType;
    // Identify the root of the given AST
    switch (expr->opType) {
        case hsql::OperatorType::kOpPlus:
            eType = ExpressionType::Addition;
            break;
        case hsql::OperatorType::kOpMinus:
            eType = ExpressionType::Subtraction;
            break;
        case hsql::OperatorType::kOpAsterisk:
            eType = ExpressionType::Multiplication;
            break;
        case hsql::OperatorType::kOpSlash:
            eType = ExpressionType::Division;
            break;
        case hsql::OperatorType::kOpEquals:
            eType = ExpressionType::Equal;
            break;
        case hsql::OperatorType::kOpNotEquals:
            eType = ExpressionType::NotEqual;
            break;
        case hsql::OperatorType::kOpLess:
            eType = ExpressionType::LessThan;
            break;
        case hsql::OperatorType::kOpLessEq:
            eType = ExpressionType::LessThanEqual;
            break;
        case hsql::OperatorType::kOpGreater:
            eType = ExpressionType::GreaterThan;
            break;
        case hsql::OperatorType::kOpGreaterEq:
            eType = ExpressionType::GreaterThanEqual;
            break;
        case hsql::OperatorType::kOpAnd:
            eType = ExpressionType::And;
            break;
        case hsql::OperatorType::kOpOr:
            eType = ExpressionType::Or;
            break;
        case hsql::OperatorType::kOpIn: // Expression is a semi-join
            eType = ExpressionType::In;
            break;
        case hsql::OperatorType::kOpExists: // Expression is an EXISTS
            eType = ExpressionType::Exists;
            break;
        default:
            eType = ExpressionType::Unknown;
    }

    if (eType == ExpressionType::Unknown) { // Expression is not an operation
        switch (expr->type) {
            case hsql::kExprColumnRef: // Expression is a column name
            {
                std::string table = expr->table == nullptr ? "" : expr->table;
                return std::make_shared<AttributeExpr>(AttributeExpr(expr->name,
                                                                  expr->table,
                                                                  expr->alias));
            }
            case hsql::kExprLiteralString: // Expression is a string
                return std::make_shared<ConstantExpr>(ConstantExpr(expr->name));
            case hsql::kExprLiteralInt: // Expression is an integer
                return std::make_shared<ConstantExpr>(ConstantExpr(expr->ival));
            case hsql::kExprLiteralFloat: // Expression is a float
                return std::make_shared<ConstantExpr>(ConstantExpr(expr->fval));
            case hsql::ExprType::kExprStar: // Expression is a '*' (all)
                return std::make_shared<Expression>(Expression(ExpressionType::Star));
            case hsql::kExprFunctionRef: // Expression is a function
            {
                auto func = std::make_shared<FunctionExpr>(FunctionExpr(expr->name,
                                                                      expr->alias,
                                                                      expr->distinct));
                // Collect input of the aggregation function
                func->expressions.push_back(processExpression(expr->exprList->at(0), query));
                return func;
            }
            case hsql::kExprSelect: // Expression is a sub-query with alias
            {
                assert(expr->alias != nullptr);
                query.subTrees.push_back(expr->select); // Update list of sub-queries
                return std::make_shared<RelationExpr>(RelationExpr(expr->alias,
                                                                   expr->alias,
                                                                   true));
            }
            default:
                leaderLogging("[ERROR] Unsupported SQL expression.\n");
                exit(-1);
        }
    }

    // Initialize expression at the current level
    auto expression = std::make_shared<Expression>(Expression(eType));

    // Collect sub-expressions (operands)
    if (expr->expr != nullptr) { // First operand
        expression->addExpression(processExpression(expr->expr, query));
    }
    if (expr->expr2 != nullptr) { // First operand
        expression->addExpression(processExpression(expr->expr2, query));
    }
    if (expr->exprList != nullptr) { // N-ary operation or sub-query with alias
        for (auto exp : expr->exprList[0]) {
            auto child = processExpression(exp, query);
            expression->addExpression(child);
        }
    }
    if (expr->select != nullptr) { // Sub-query without alias
        int s = query.subTrees.size();
        // Anonymous sub-queries are assigned an id of the form "#TMP_"
        std::string name = "#TMP_"+std::to_string(s);
        query.subTrees.push_back(expr->select); // Update list of sub-queries
        auto child = std::make_shared<RelationExpr>(RelationExpr(name.c_str(),
                                                                                        nullptr,
                                                                                        true));
        expression->addExpression(child);
    }

    // Make sure the expression has been parsed correctly
    assert(eType!=ExpressionType::Unknown);

    return expression;
}

// Processes ORDER BY clause
static void processWithClause(const SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc) {
    std::cout << "COMPILES" << std::endl;
}

#endif //SECRECY_PARSER_UTILS_H
