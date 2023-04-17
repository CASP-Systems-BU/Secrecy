#ifndef SECRECY_PLANGEN_H
#define SECRECY_PLANGEN_H

#include "plan.h"
#include "rule.h"
#include "queryGraph.h"
#include "dpTable.h"
#include "database.h"
#include "bitset.h"
#include "costs.h"

// A logical plan generator 
class PlanGen {
private:
    // The DP table
    std::vector<DPTable> dpts;
    // The database
    Database *db;
    // The current query
    std::shared_ptr<QueryGraph> qGraph;

    // Apply scan
    std::shared_ptr<Problem> applyScan(unsigned qId, std::shared_ptr<RelationExpr> rel);
    // Apply project
    void applyProject(unsigned qId, const std::vector<std::shared_ptr<Expression>> &projections,
                      std::shared_ptr<Plan> plan);
    // Apply filter
    std::shared_ptr<Plan> applyFilter(std::shared_ptr<FilterExpr> filter, std::shared_ptr<Plan> plan);
    // Apply join
    std::shared_ptr<Plan> applyJoin(std::shared_ptr<JoinExpr> join, std::shared_ptr<Plan> left,
                                    std::shared_ptr<Plan> right);
    // Apply join as a filter
    std::shared_ptr<Plan> applyJoin(std::shared_ptr<JoinExpr> join, std::shared_ptr<Plan> plan);
    // Apply semi-join
    std::shared_ptr<Plan> applySemiJoin(std::shared_ptr<JoinExpr> join, std::shared_ptr<Plan> left,
                                        std::shared_ptr<Plan> right);
    // Apply distinct
    void applyDistinct(std::shared_ptr<QueryGraph> query, std::shared_ptr<Problem> prob);
    // Apply group-by
    void applyGroupBy(std::shared_ptr<QueryGraph> query, std::shared_ptr<Problem> prob);
    // Apply order-by
    void applyOrderBy(std::shared_ptr<QueryGraph> query, std::shared_ptr<Problem> prob);
    // Apply table function
    std::shared_ptr<Plan> applyFunction(std::shared_ptr<TableFunctionExpr> function,
                                        std::shared_ptr<Plan> plan);

    // Initialize DP table with scans
    void applyScans(unsigned qId, std::vector<std::shared_ptr<RelationExpr>> &relations);
    // Expand solutions to the given problem with applicable filters
    void applyFilters(unsigned qId, std::vector<std::shared_ptr<FilterExpr>> &filters,
                      std::shared_ptr<Problem> prob);
    // Expand solutions to the given problem with applicable joins
    void applyJoins(unsigned qId, std::vector<std::shared_ptr<JoinExpr>> &joins,
                    std::shared_ptr<Problem> prob);
    // Combine solutions to the given pair of problems with applicable joins
    void applyJoins(unsigned qId, std::vector<std::shared_ptr<JoinExpr>> &joins,
                    std::shared_ptr<Problem> prob1, std::shared_ptr<Problem> prob2);
    // Expand solutions to the given problem with applicable functions
    void applyFunctions(std::shared_ptr<QueryGraph> query, std::shared_ptr<Problem> prob);
    // Apply all possible unary operators to the plans of the given problem
    void applyUnaryOperators(unsigned qId, std::shared_ptr<QueryGraph> qGraph,
                             std::shared_ptr<Problem> prob);

    // Pick a candidate final plan
    std::shared_ptr<Plan> pickOptimalPlan(unsigned qId);
    // Apply logical optimization rules
    void applyRules(int qId);
    // Refines expressions for the code generator
    void updateAttributes(std::shared_ptr<Plan> op, std::vector<std::shared_ptr<Expression>> &exp);
    // Refines plan metadata for the code generator
    void postProcess(std::vector<std::shared_ptr<Plan>> &plans);
public:
    /***
     * Constructor
     */
    PlanGen() { };
    /***
     * Destructor
     */
    ~PlanGen() { };

   /***
    * Translate a query graph into a logical execution plan.
    * @param qGraph - The input query graph.
    * @param plans - The logical plan(s). May be more than one in case of nested queries.
    * @param db - The database the query is evaluated on.
    * @return The optimal logical plan for the given query.
    */
   std::shared_ptr<Plan> generatePlan(std::shared_ptr<QueryGraph> qGraph, std::vector<std::shared_ptr<Plan>> &plans,
                                      Database &db);
};

#endif  // SECRECY_PLANGEN_H