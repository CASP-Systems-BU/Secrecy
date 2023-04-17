#ifndef SECRECY_PREPROCESSOR_H
#define SECRECY_PREPROCESSOR_H

#include "plan.h"
#include "queryGraph.h"
#include "database.h"
#include "bitset.h"
#include "costs.h"
#include <queue>
#include "sql_expression.h"

// The prefix name for anonymous sub-queries
const std::string TMP_PREFIX = "#Q_";

// A query preprocessor
class Preprocessor {
private:
    // Number of (sub-)queries processed
    int qc = 0;
    // The queue of queries to process
    std::queue<std::pair<Parser::SQLQuery, std::shared_ptr<QueryGraph>> > queue;
    // The query index
    std::map<int, std::shared_ptr<QueryGraph>> qIndex;
    // The sub-query index
    std::map<const hsql::SelectStatement*, std::vector<std::shared_ptr<RelationExpr>>> sIndex;
    // The aliases-to-id mapping
    std::map<std::string, int> a2id;

    // Checks if state is empty
    bool stateIsEmpty() const;
    // Sets the name of the temporary relation that represents the output of the given sub-query
    void setAlias(Parser::SQLQuery &subQuery);
    // Retrieves the graph of the query with the given `name`
    std::shared_ptr<QueryGraph> resolveQuery(std::shared_ptr<RelationExpr> rel) const;
    // Resolves the result schema of the given query
    void resolveSchema(std::shared_ptr<QueryGraph> query, Database &db,
                       std::shared_ptr<RelationExpr> relation=nullptr) const;
    // Checks if the given attribute exists in the given (base or temporary) relation
    bool validate(std::shared_ptr<AttributeExpr> att, std::shared_ptr<RelationExpr> rel, Database &db) const;
    // Checks if a relation with name `rName` exists in the database
    void validate(std::string rName, Database &db) const;
    // Checks if all base relations in the query exist in the database
    void validateRelations(std::shared_ptr<QueryGraph> query, Database &db) const;
    // Check if the given attribute expression is valid w.r.t. the given query and database
    void validate(std::shared_ptr<AttributeExpr> a, std::shared_ptr<QueryGraph> query, Database &db,
                  bool sortingKey=false) const;
    // Validates all attributes in the given expression
    void validate(std::shared_ptr<Expression>& expr, std::shared_ptr<QueryGraph> &query, Database &db,
                  bool sortingKey=false) const;
    // Validates attribute names and updates their relation information accordingly
    void validateAttributes(std::shared_ptr<QueryGraph> query, Database &db) const;
    // Validates IN operators in the given expression
    void validateSemiJoins(std::shared_ptr<Expression> expr, std::shared_ptr<QueryGraph> query,
                           Database &db) const;
    // Validates IN operators in the query
    void validateSemiJoins(std::shared_ptr<QueryGraph> query, Database &db) const;

    // Processes a parsed SQL expression
    std::shared_ptr<Expression> processExpression(const hsql::Expr *expr, Parser::SQLQuery& query,
                                                  std::shared_ptr<Expression> parent=nullptr);
    // Processes a reference to base relation
    std::shared_ptr<RelationExpr> processTableRef(const hsql::TableRef* tableRef);
    // Processes a join reference
    std::shared_ptr<JoinExpr> processJoin(const hsql::JoinDefinition* joinRef, Parser::SQLQuery &query);
    // Processes SELECT clause
    void processSelectClause(Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc);
    // Processes FROM clause
    void processFromClause(Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc);
    // Processes WHERE clause
    void processWhereClause(Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc);
    // Processes GROUP BY clause
    void processGroupByClause(Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc);
    // Processes ORDER BY clause
    void processOrderByClause(const Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc);
    // Processes WITH clause
    void processWithClause(const Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc);
    // Processes queries in the queue
    void processQueue();
    // Collects information about the given query
    std::shared_ptr<QueryGraph> processQuery(Parser::SQLQuery& query, int id=0,
                                             std::shared_ptr<QueryGraph> parent=nullptr);
    // Collects information about the given query and its sub-queries (if any)
    std::shared_ptr<QueryGraph> collectInfo(Parser::SQLQuery& query);

public:
    // Constructor
    Preprocessor();
    // Destructor
    ~Preprocessor();

    /***
     * Applies basic query validation, e.g., checks for non-existent or ambiguous columns, invalid operations, etc.
     * @param query - The query to validate.
     * @param db - The database the query is applied to.
     */
    void validate(std::shared_ptr<QueryGraph> query, Database &db);

    /***
     * Applies semantic analysis and transforms the given query into an equivalent one, e.g., by removing duplicate
     * predicates, rewriting nested queries, etc.
     * @param query - The query to preprocess.
     * @param db - The database the query is applied to.
     */
    void preprocess(std::shared_ptr<QueryGraph> query, Database &db);

    /***
     * Transforms a parsed SQL query into a query graph.
     * @param query - The parsed SQL query.
     * @param db - The database the query is applied to.
     * @param qGraph - The query graph to populate.
     */
    std::shared_ptr<QueryGraph> transform(Parser::SQLQuery &query, Database &db);
};

#endif  // SECRECY_PREPROCESSOR_H