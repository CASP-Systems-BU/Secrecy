#include "../../include/planner/dpTable.h"

void Problem::printSolutions() {
    if (get_rank()==0) {
        for (auto it=plans; it!=nullptr; it=it->next) {
            it->print(0);
        }
    }
}

std::shared_ptr<Problem> DPTable::finalProblem() {
    unsigned stage = stages() - 1;
    unsigned entry = entries(stage) - 1;
    auto pr = dptable[stage][entry];
    assert(pr != nullptr);
    return pr;
}

void DPTable::printSolutions() const {
    if (get_rank()==0) {
        unsigned cnt=0;
        for (int i=0; i<this->stages(); i++) { // For all stages
            for (int j=0; j<this->entries(i); j++) { // For all entries per stage
                for (auto pr=this->dptable[i][j]; pr!=nullptr; pr=pr->next) { // For all problems
                    for (auto pl=pr->plans; pl!=nullptr; pl=pl->next) { // For all plans per problem
                        std::cout << "Plan " << cnt++ << " (Cost: " << pl->absoluteCost() << "):" << std::endl;
                        pl->print(0);
                        std::cout << "========================" << std::endl;
                    }
                }
            }
        }
    }
}

void DPTable::add(std::shared_ptr<Problem> prob, unsigned pos) {
    assert((prob != nullptr) && (this->dptable.size() > pos));
    if (this->dptable[pos][0] != nullptr)
        prob->next = this->dptable[pos][0];
    this->dptable[pos][0] = prob;
}

// Adds a plan to the set of possible solutions to the problem
bool DPTable::addPlan(std::shared_ptr<Problem> problem, std::shared_ptr<Plan> plan) {
    int i=0;
    std::shared_ptr<Plan> last = nullptr;
    for (std::shared_ptr<Plan> pl = problem->plans, next; pl!=nullptr; pl=pl->next) {
        if (pl->dominates(plan)) { // Candidate plan is dominated by existing plan
            return false;
        }
        if (plan->dominates(pl)) { // Candidate plan dominates existing plan
            if (last != nullptr) {
                last->next = pl->next;
            }
            else {
                problem->plans = pl->next;
            }
        }
        else {
            last = pl;
        }
        i++;
    }
    if (last) { // Append plan to the list of possible solutions
        last->next = plan;
    }
    else { // First plan
        plan->next = problem->plans;
        problem->plans = plan;
    }
    return true;
}

void DPTable::add(std::shared_ptr<Plan> plan, BitSet& relations) {
    assert(plan != nullptr);
    unsigned stage = plan->ops, entry = plan->nodes - 1;
    assert((stages() > stage) && (entries(stage) > entry));
    for (auto pr=dptable[stage][entry]; pr!=nullptr; pr=pr->next) {
        if (pr->relations==relations) {
            this->addPlan(pr, plan);
            return;
        }
    }
    auto prob = std::make_shared<Problem>(Problem());
    prob->relations.set(relations);
    if (this->dptable[stage][plan->nodes-1] != nullptr)
        prob->next = this->dptable[stage][0];
    prob->plans = plan;
    this->dptable[stage][entry] = prob;
}

std::shared_ptr<Problem> DPTable::at(unsigned stage, unsigned entry) {
    return this->dptable[stage][entry];
}

unsigned DPTable::stages() const {
    return this->dptable.size();
}

unsigned DPTable::entries(unsigned stage) const {
    return this->dptable[stage].size();
}

void DPTable::allocate(unsigned qId, int stages, int entries) {
    this->dptable.resize(stages);
    for (int i=0; i<this->dptable.size(); i++)
        this->dptable[i].resize(std::min(i+1, entries));
    // Print info
    std::string msg = "[INFO] Initialized container with " + std::to_string(this->dptable.size()) + " slots "
                      + " for query with id " + std::to_string(qId) + "\n";
    leaderLogging(msg);
    leaderLogging("Entries per slot: ");
    msg.clear();
    for (int i=0; i<this->dptable.size(); i++)
        msg += std::to_string(this->dptable[i].size()) + " ";
    leaderLogging(msg+"\n");
}