#ifndef SECRECY_QUERYGRAPH_H
#define SECRECY_QUERYGRAPH_H

#include <vector>
#include <string>
#include <set>
#include <map>
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <climits>

#include "database.h"
#include "parser/parser.h"
#include "relation.h"
#include "sql_expression.h"

/***
 * A query graph representing a SQL query.
 */
class QueryGraph : public std::enable_shared_from_this<QueryGraph>{
private:
    // The cardinalities type
    typedef double card_t;

    // Checks if the given expression is a filter and returns the respective filter expression
    std::shared_ptr<FilterExpr> isFilter(const std::shared_ptr<Expression> expr) const;
    // Checks if the given expression is a join and returns the respective join expression
    std::shared_ptr<JoinExpr> isJoin(std::shared_ptr<Expression> expr,
                                     JoinExpr::JoinType joinType=JoinExpr::Inner) const;
    // Collects information about filters and joins in the given expression
    void collectFiltersAndJoins(std::shared_ptr<Expression> expr, JoinExpr::JoinType joinType=JoinExpr::Inner);
    // Collects information about filters and joins in the query
    void collectFiltersAndJoins();
    // Collects information about global aggregations in the query
    void collectTableFunctions();
    // Update base relation cardinalities
    void setBaseCardinalities(const Database& db);
    // Sets logical operator ids
    void setOpIds();
    // Prints query info
    void printInfo() const;
    // Prints query hierarchy
    void printHierarchy(int indent=0) const;

public:
    /***
     * Query ID
     */
    int id;
    /***
     * Max number of rows in the result
     */
    unsigned limit;
    /***
     * SELECT DISTINCT flag
     */
    bool distinct;
    /***
     * True if query has been validated, False otherwise
     */
    bool validated = false;
    /***
     * Parent query
     */
    std::shared_ptr<QueryGraph> parent = nullptr;
    /***
     * Sub-queries
     */
    std::vector<std::shared_ptr<QueryGraph>> children;
    /***
     * Sibling query (e.g., in case of WITH or UNION)
     */
    std::shared_ptr<QueryGraph> sibling = nullptr;
    /***
     * The list of attributes and expressions in the SELECT clause
     */
    std::vector<std::shared_ptr<Expression>> projections;
    /***
     * The list of relations involved in the query
     */
    std::vector<std::shared_ptr<RelationExpr>> relations;
    /***
     * The expression tree in the WHERE clause
     */
    std::shared_ptr<Expression> whereClause = nullptr;
    /***
     * The expression tree in the FROM clause
     */
    std::shared_ptr<Expression> fromClause = nullptr;
    /***
     * The type of join in the FROM clause
     */
    JoinExpr::JoinType jType;
    /***
     * The list of filters in the query
     */
    std::vector<std::shared_ptr<FilterExpr>> filters;
    /***
     * The list of joins in the query
     */
    std::vector<std::shared_ptr<JoinExpr>> joins;
    /***
     * Group-by clause
     */
    std::shared_ptr<GroupByExpr> groupBy = nullptr;
    /***
     * Order-by clause
     */
    std::shared_ptr<OrderByExpr> orderBy = nullptr;
    /***
     * The table functions (e.g., global aggregations) in the query
     */
    std::vector<std::shared_ptr<TableFunctionExpr>> tFunctions;
    /***
     * The temporary relation names `this` (sub-)query corresponds to
     */
    std::vector<std::shared_ptr<RelationExpr>> relAlias;

    /***
     * Constructor
     * @param id - The id of the query
     * @param parent - A pointer to the parent query (if any)
     */
    QueryGraph(int id, std::shared_ptr<QueryGraph> parent=nullptr) : id(id), parent(parent) { };
    /***
     * Destructor
     */
    ~QueryGraph() { };

    /***
     * @param a - An attribute expression.
     * @return True if the query output includes the given attribute, False otherwise.
     */
    bool outputs(std::shared_ptr<AttributeExpr> a);
    /***
     * @return True if the query contains a '*' (all)
     */
    bool containsStar() const;
    /***
     * @return True if `this` query appears in parent's WHERE clause, False otherwise
     */
    bool appearsInWhereClause() const;
    /***
     * Constructs a mapping from query IDs to (sub-)queries. Queries are numbered in a BFS fashion.
     * @param m - The map to populate.
     */
    void collectQueries(std::map<int, std::shared_ptr<QueryGraph>> &m);
    /***
    * Initializes the query graph
    */
    void constructEdges(const Database& relations);


    /***
     * Retrieve all expressions that correspond to the given relation
     * @param rel_id - The id of the relation whose expressions we want to retrieve.
     * @param v - A vector to populate.
     * @param where - If True, the method will also search in the WHERE clause of the query for any temporary relations
     */
    void getRelations(const std::string &rel_id, std::vector<std::shared_ptr<RelationExpr> > &v, bool where=false);

    /***
     *
     * @param rel_id
     * @param where
     * @return
     */
    std::shared_ptr<RelationExpr> getRelation(const std::string &rel_id, bool where=false) const;

    /***
     * Returns the (integer) id of the relation with the given name or alias
     * @param id
     * @return
     */
    int getRelId(const std::string &id) const;


    /***
     * Checks if the given relation exists in `this` query.
     * @param rel_name - The name of the relation to search for.
     * @return True if the relation exists in the query, False otherwise.
     */
    bool contains(const std::string &rel_name) const;

    /***
     * @return True if `this` query includes a global or per-group aggregation, False otherwise.
     */
    bool hasAggregation() const;
    /***
     * @return The number of logical operators in `this` query.
     */
    unsigned numOp() const;
    /***
     * @return The id of the distinct operator.
     */
    unsigned getDistinctId() const;
    /***
     * @return True if `this` query has a parent.
     */
    bool isSubQuery() const { return !relAlias.empty(); };
    /***
     * Prints complete query information
     */
    void print() const;

};
#endif // SECRECY_QUERYGRAPH_H