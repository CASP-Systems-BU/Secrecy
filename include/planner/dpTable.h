#ifndef SECRECY_DPTABLE_H
#define SECRECY_DPTABLE_H

#include <algorithm>

#include "rule.h"
#include "plan.h"
#include "database.h"
#include <iterator>
#include <cstddef>

// A problem
struct Problem {
    // The next problem in the DP table
    std::shared_ptr<Problem> next;
    // The known solutions to the problem
    std::shared_ptr<Plan> plans;
    // The relations involved in the problem
    BitSet relations;
    // Prints known solutions
    void printSolutions();
};

class DPTable {
public:
    // The DP table
    std::vector<std::vector<std::shared_ptr<Problem>>> dptable;

    void allocate(unsigned qId, int stages, int entries);

    void add(std::shared_ptr<Problem> prob, unsigned pos);

    void add(std::shared_ptr<Plan> plan, BitSet& relations);

    /***
     * Retrieves problem at the given `stage` and `entry`
     * @param stage - The stage.
     * @param entry - The entry.
     * @return The first problem at the given `entry`.
     */
    std::shared_ptr<Problem> at(unsigned stage, unsigned entry);

    /***
     * Retrieves the total number of stages in the DP table. A stage is identified by the number of operators included
     * in its plans, i.e., stage 0 includes scans, stage 1 includes plans with one operator, etc.
     * @return The total number of stages in the table.
     */
    unsigned stages() const;
    /***
     * Retrieves the total number of entries for the given stage of the DP table. An entry is identified by the number
     * of relations included in its plans, i.e., entry 0 includes all plans with 1 relation, entry 1 includes all plans
     * with 2 relations, etc.
     * @param stage - The stage whose entries we want to retrieve.
     * @return The number of entries for the given stage.
     */
    unsigned entries(unsigned stage) const;

    /***
     * Adds `plan` to the set of possible solutions to the `problem`.
     */
    bool addPlan(std::shared_ptr<Problem> problem, std::shared_ptr<Plan> plan);

    /***
     * @return The problem from the last stage of the DP table.
     */
    std::shared_ptr<Problem> finalProblem();

    /***
     * Prints all (sub-)plans in the DP table to stdout.
     */
    void printSolutions() const;
};

#endif //SECRECY_DPTABLE_H
