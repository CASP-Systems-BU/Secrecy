#ifndef SECRECY_PLAN_H
#define SECRECY_PLAN_H

#include <iostream>
#include <stack>

#include "bitset.h"
#include "costs.h"
#include "database.h"
#include "operators/expression.h"


class SQLExpression;

/***
 * A logical plan
 */
class Plan : public std::enable_shared_from_this<Plan> {
private:
    // Collects operators in `this` plan
    void collectOps(std::vector<std::shared_ptr<Plan>> &ops);
public:
    /***
     * The cardinalities type
     */
    typedef double card_t;
    /***
     * The absolute cost type
     */
    typedef double cost_t;
    /***
     * Operator types
     */
    enum Op { Scan,
            NestedLoopJoin,
            SemiJoin,
            Filter,
            JoinFilter,
            Project,
            Union,
            TableFunction,
            GroupBy,
            OddEvenMerge,
            Distinct,
            OrderBy,
            GroupJoin,
            RowFunction,
            AdjacentEq,
            AdjacentGeq,
            Composition,
            Mask};

    /***
     * Constructor
     */
    Plan();
    /***
     * Copy constructor
     */
    Plan(const std::shared_ptr<const Plan> other);
    /***
     * Destructor
     */
    ~Plan();

    /***
     * The root operator type
     */
    Op op;
    /***
     * The sub-plans
     */
    std::shared_ptr<Plan> left=nullptr, right=nullptr;
    /***
     * The output cardinality
     */
    card_t cardinality;
    /***
     * The exact plan cost
     */
    CostModel::Cost planCost;
    /***
     * The exact predicate cost (for filters and joins)
     */
    CostModel::Cost predCost;
    /***
     * The ordering, i.e., the id of the attribute whose values are in order (if any)
     */
    unsigned ordering = UINT32_MAX;
    /***
     * The next plan for the same problem
     */
    std::shared_ptr<Plan> next;
    /***
     *  The id of the base relation (for scans)
     */
    unsigned node = UINT32_MAX;
    /***
     *  The id of the root operator
     */
    unsigned opId = UINT32_MAX;
    /***
     * Total number of relations in the plan
     */
    unsigned nodes = 0;
    /***
     * Total number of physical operators in the plan
     */
    unsigned ops = 0;
    /***
     * A mask of operators in the plan
     */
    BitSet operators;
    /***
     * The SQL expression the root operator corresponds to
     */
    std::shared_ptr<SQLExpression> exp;
    /***
     * Projected attributes
     */
    std::vector<std::shared_ptr<Expression>> attributes;

    /***
     *  Checks if `this` plan costs less than the `other` plan
     */
    bool dominates(const std::shared_ptr<Plan> other) const;
    /***
     * Checks if the plan contains the given expression.
     * @param op - The expression to search for.
     * @return True if `op` is included in `this` plan, False otherwise.
     */
    bool contains(std::shared_ptr<SQLExpression> op) const;
    /***
     * Returns the absolute plan cost
     */
    unsigned absoluteCost() const;
    /***
     *  Returns true if the root operator is a join
     */
    bool rootIsJoin() const;
    /***
     * Retrieves the base relations involved in the plan
     */
    void baseRelations(std::vector<std::shared_ptr<RelationExpr>> &relations) const;
    /***
     * Checks if the given attributes belong to one input relation
     * @param attributes - The set of attributes to search for.
     * @return 0 if attributes come from left, 1 if they come from right, -1 if they come from multiple relations
     */
    int singleOrigin(const std::vector<std::shared_ptr<Expression>> &attributes) const;
    /***
     * Checks if the root operator of `this` plan is a self join.
     * @param adjacent - A flag to check whether the self join is on adjacent rows (True) or not (False).
     * @return True if the self-join pattern is found, False otherwise.
     */
    bool isSelfJoin(bool adjacent=false) const;
    /***
     * Checks if the root operator of `this` plan is a join and retrieves its predicates.
     * @param predicates - The vector of predicates to populate.
     * @return True if the root of `this` plan is a join, False otherwise.
     */
    bool isJoin(std::vector<std::shared_ptr<SQLExpression>> &predicates) const;
    /***
     * Checks if the root operator of `this` plan is a distinct
     */
    bool isDistinct() const;
    /***
     * Retrieves the scan operators in the plan.
     * @param scans - The vector of plans to populate.
     */
    void collectScans(std::vector<std::shared_ptr<Plan>> &scans);

    void getJoinInputs(std::vector<std::shared_ptr<Plan>> &plans);
    /***
     * Retrieves all plan operators in topological order (i.e., scans first, operators above scans, ..., root)
     * @param operators - The vector of operators (sub-plans) to populate.
     */
    void collectOperators(std::vector<std::shared_ptr<Plan>> &operators);
    /***
     * @return True, if `this` plan requires a composition operator, False otherwise.
     */
    bool requiresComposition() const { return (op != Scan && op != OrderBy); };
    /***
     * @return The derived attribute that is appended by the root operator of `this plan`.
     */
    std::shared_ptr<Expression> getDerivedAttribute();
    /***
     * Retrieves the indices of the key attributes.
     * @param key_ops - The vector of indices to populate.
     */
    void getKeyPositions(std::vector<unsigned> &key_ops) const;
    /***
     * @return True if this plan returns a derived arithmetic attribute
     */
    bool hasAAttribute() const;
    /***
     *  Prints the plan to stdout
     */
    void print(unsigned indent=0, bool pCost=false, bool rootOnly=false, bool schema=false) const;
};

#endif // SECRECY_PLAN_H