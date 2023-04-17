#ifndef SECRECY_CODEGEN_H
#define SECRECY_CODEGEN_H

#include "plan.h"
#include "../core/relational.h"
#include "database.h"
#include "planGen.h"
#include "../planner/operators/operators.h"

// A physical plan generator
class CodeGen {
private:
    // The database
    Database *db;

    // Retrieves index of input operator
    unsigned getInputOp(std::shared_ptr<Plan> inputOp, const std::vector<std::shared_ptr<Plan>> &operators);
    // Retrieves the relation index
    unsigned getRelationInd(const std::vector<std::shared_ptr<Relation>> &relations, std::string name);
    // Generates physical operators
    unsigned generateOperators(std::vector<std::shared_ptr<Plan>> &operators,
                               const std::vector<std::shared_ptr<Relation>> &relations,
                               std::vector<std::shared_ptr<Operator>> &plan_ops);
public:
    /***
     * Constructor
     */
    CodeGen() {};
    /***
     * Destructor
     */
    ~CodeGen() {};

    /***
      * Generates the physical plan for the given logical `plan` and database `db`
      * @param plan - The logical plan.
      * @param db - The database the query is applied to.
      * @return The root operator of the physical plan.
      */
    std::shared_ptr<Operator> generatePhysicalPlan(std::vector<std::shared_ptr<Plan>> &plan, Database *db);
};

#endif //SECRECY_CODEGEN_H
