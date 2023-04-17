#include "../../include/planner/planGen.h"

std::shared_ptr<Problem> PlanGen::applyScan(unsigned qId, std::shared_ptr<RelationExpr> rel) {
    assert(!rel->schema.empty());
    // Create a new problem
    std::shared_ptr<Problem> prob = std::make_shared<Problem>(Problem());
    prob->next = nullptr;
    prob->plans = nullptr;
    prob->relations = BitSet();
    prob->relations.set(rel->id);
    // Initialize a new plan
    auto scan = std::make_shared<Plan>(Plan());
    scan->op = Plan::Scan;
    auto scExpr = std::make_shared<ScanExpr>(ScanExpr());
    scExpr->rel = rel;
    scan->exp = scExpr;
    scan->node = rel->id;
    scan->nodes = 1;
    scan->cardinality = rel->cardinality;
    for (auto a : rel->schema) {
        if (a->relation != rel->name) { // Check for aliases
            if (auto atr = std::dynamic_pointer_cast<AttributeExpr>(a)) {
                auto att = std::make_shared<AttributeExpr>(AttributeExpr(atr));
                rel->alias.empty() ? att->relation = rel->name : att->relation = rel->alias;
                scan->attributes.push_back(att);
                continue;
            }
        }
        scan->attributes.push_back(a);
    }


    scan->planCost = CostModel::scanCost(scan->cardinality);
    // Add plan to the set of possible solutions
    this->dpts[qId].addPlan(prob, scan);
    return prob;
}

std::shared_ptr<Plan> PlanGen::applyFilter(std::shared_ptr<FilterExpr> filter, std::shared_ptr<Plan> plan) {
    // The new plan after the filter application
    auto np = std::make_shared<Plan>(Plan());
    np->op = Plan::Filter;
    np->exp = filter;
    np->opId = filter->id;
    np->left = plan;
    np->cardinality = plan->cardinality;
    np->ordering = plan->ordering;
    np->planCost = plan->planCost + filter->cost(plan->cardinality);;
    // TODO: Composition
    for (auto a : plan->attributes) {
        if (!a->hasType(Expression::DerivedAttributeRef))
            np->attributes.push_back(a);
    }
    np->attributes.push_back(std::make_shared<DerivedAttributeExpr>("S_BIT", filter->id));
    np->nodes = plan->nodes;
    np->node = plan->node;
    np->ops = plan->ops + 1;
    np->operators = plan->operators;
    np->operators.set(np->opId);
    return np;
}

void PlanGen::applyProject(unsigned qId, const std::vector<std::shared_ptr<Expression>> &projections,
                           std::shared_ptr<Plan> plan) {
    std::vector<std::shared_ptr<Expression>> new_proj;
    for (auto a : plan->attributes) {
        if (a->hasType(Expression::DerivedAttributeRef)) {
            new_proj.push_back(a);
            continue;
        }
        if (auto att = std::dynamic_pointer_cast<AttributeExpr>(a)) {
            if (att->name.find("SEL") != std::string::npos) {
                new_proj.push_back(a);
                continue;
            }
        }
        for (auto p : projections) {
            if (a->sameAs(p)) {
                new_proj.push_back(p);
                break;
            }
        }
    }
    // Update output schema
    plan->attributes = new_proj;
}

std::shared_ptr<Plan> PlanGen::applyJoin(std::shared_ptr<JoinExpr> join, std::shared_ptr<Plan> left,
                                std::shared_ptr<Plan> right) {
    // The new plan after the join application
    auto np = std::make_shared<Plan>(Plan());
    np->op = Plan::NestedLoopJoin;
    np->opId = join->id;
    np->exp = join;
    np->left = left;
    np->right = right;
    np->cardinality = left->cardinality * right->cardinality;
    np->ordering = left->ordering;
    // TODO: Composition
    np->predCost = join->cost();
    np->planCost = left->planCost + right->planCost + CostModel::joinCost(left->cardinality,
                                                                          right->cardinality,
                                                                          np->predCost);
    np->nodes = left->nodes + right->nodes;
    np->node = UINT32_MAX;
    np->ops = left->ops + right->ops + 1;
    np->operators = left->operators.unionWith(right->operators);
    np->operators.set(np->opId);
    np->attributes = left->attributes;
    np->attributes.insert(np->attributes.end(), right->attributes.begin(), right->attributes.end());
    np->attributes.push_back(std::make_shared<DerivedAttributeExpr>("J_BIT", join->id));
    return np;
}

std::shared_ptr<Plan> PlanGen::applySemiJoin(std::shared_ptr<JoinExpr> join, std::shared_ptr<Plan> left,
                                             std::shared_ptr<Plan> right) {
    // The new plan after the semi-join application
    auto np = std::make_shared<Plan>(Plan());
    np->op = Plan::NestedLoopJoin;
    np->opId = join->id;
    np->exp = join;
    np->left = left;
    np->right = right;
    np->cardinality = left->cardinality;
    np->ordering = left->ordering;
    // TODO: Composition
    np->attributes = left->attributes;
    np->attributes.push_back(std::make_shared<DerivedAttributeExpr>("SJ_BIT", join->id));
    np->predCost = join->cost();
    np->planCost = left->planCost + right->planCost + CostModel::semiJoinCost(left->cardinality,
                                                                              right->cardinality,
                                                                              np->predCost);
    np->nodes = left->nodes + right->nodes;
    np->node = UINT32_MAX;
    np->ops = left->ops + right->ops + 1;
    np->operators = left->operators.unionWith(right->operators);
    np->operators.set(np->opId);
    return np;
}

std::shared_ptr<Plan> PlanGen::applyJoin(std::shared_ptr<JoinExpr> join, std::shared_ptr<Plan> plan) {
    // The new plan after the join application
    auto np = std::make_shared<Plan>(Plan());
    np->op = Plan::JoinFilter; // Join is reduced into a filter
    np->opId = join->id;
    np->exp = join;
    np->left = plan;
    np->cardinality = plan->cardinality;
    np->ordering = plan->ordering;
    // TODO: Composition
    for(auto a : plan->attributes)
        if (!a->hasType(Expression::DerivedAttributeRef))
            np->attributes.push_back(a);
    np->attributes.push_back(std::make_shared<DerivedAttributeExpr>("J_BIT", join->id));
    np->predCost = join->cost();
    np->planCost = plan->planCost + join->cost(np->cardinality);
    np->nodes = plan->nodes;
    np->node = UINT32_MAX;
    np->ops = plan->ops + 1;
    np->operators = plan->operators;
    np->operators.set(np->opId);
    return np;
}

void PlanGen::applyDistinct(std::shared_ptr<QueryGraph> query, std::shared_ptr<Problem> prob) {
    if (!query->distinct)
        return;
    unsigned n = query->filters.size() + query->joins.size() + (query->groupBy != nullptr);
    for (auto pl=prob->plans; pl!=nullptr; pl=pl->next) {
        if (pl->ops == n) {
            // The new plan after the distinct application
            auto np = std::make_shared<Plan>(Plan());
            np->op = Plan::Distinct;
            np->opId = query->getDistinctId();
            // Create distinct expression
            auto dExpr = std::make_shared<DistinctExpr>(DistinctExpr(np->opId));
            dExpr->keys = query->projections;
            np->exp = dExpr;
            np->left = pl;
            np->cardinality = pl->cardinality;
            np->ordering = pl->ordering;
            // TODO: Composition
            np->attributes = dExpr->keys;
            np->attributes.push_back(std::make_shared<DerivedAttributeExpr>("D_BIT", np->opId));
            auto c = CostModel::distinctCost(np->cardinality, pl->requiresComposition());
            np->planCost = pl->planCost + c;
            np->nodes = pl->nodes;
            np->node = UINT32_MAX;
            np->ops = pl->ops + 1;
            np->operators = pl->operators;
            np->operators.set(np->opId);
            dpts[query->id].add(np, prob->relations);
        }
    }
}

void PlanGen::applyGroupBy(std::shared_ptr<QueryGraph> query, std::shared_ptr<Problem> prob) {
    if (!query->groupBy)
        return;
    unsigned n = query->filters.size() + query->joins.size();
    for (auto pl=prob->plans; pl!=nullptr; pl=pl->next) {
        if (pl->ops == n) {
            // The new plan after the distinct application
            auto np = std::make_shared<Plan>(Plan());
            np->op = Plan::GroupBy;
            np->exp = query->groupBy;
            np->opId = query->groupBy->id;
            np->left = pl;
            np->cardinality = pl->cardinality;
            np->ordering = pl->ordering;
            np->attributes = query->projections;
            // TODO: Composition
            auto aggCost = CostModel::Cost();
            auto c = CostModel::groupByCost(np->cardinality, aggCost);
            np->planCost = pl->planCost + c;
            np->nodes = pl->nodes;
            np->node = UINT32_MAX;
            np->ops = pl->ops + 1;
            np->operators = pl->operators;
            np->operators.set(np->opId);
            dpts[query->id].add(np, prob->relations);
            np->exp = query->groupBy;
        }
    }
}

void PlanGen::applyOrderBy(std::shared_ptr<QueryGraph> query, std::shared_ptr<Problem> prob) {
    if (!query->orderBy)
        return;
    unsigned n = query->numOp() - 1;
    for (auto pl=prob->plans; pl!=nullptr; pl=pl->next) {
        if (pl->ops == n) {
            // The new plan after the distinct application
            auto np = std::make_shared<Plan>(Plan());
            np->op = Plan::OrderBy;
            np->exp = query->orderBy;
            np->opId = query->orderBy->id;
            np->left = pl;
            np->cardinality = pl->cardinality;
            np->ordering = pl->ordering;
            np->attributes = pl->attributes;
            // TODO: Composition
            auto c = CostModel::sortCost(np->cardinality);
            np->planCost = pl->planCost + c;
            np->nodes = pl->nodes;
            np->node = UINT32_MAX;
            np->ops = pl->ops + 1;
            np->operators = pl->operators;
            np->operators.set(np->opId);
            dpts[query->id].add(np, prob->relations);
        }
    }
}

std::shared_ptr<Plan> PlanGen::applyFunction(std::shared_ptr<TableFunctionExpr> function,
                                             std::shared_ptr<Plan> plan) {
    // The new plan after the filter application
    auto np = std::make_shared<Plan>(Plan());
    np->op = Plan::TableFunction;
    np->opId = function->id;
    np->exp = function;
    np->left = plan;
    np->cardinality = plan->cardinality;
    np->attributes.push_back(function->exp);
    np->ordering = plan->ordering;
    np->planCost = plan->planCost + function->cost(plan->cardinality);;
    np->nodes = plan->nodes;
    np->node = plan->node;
    np->ops = plan->ops + 1;
    np->operators = plan->operators;
    np->operators.set(np->opId);
    return np;
}

void PlanGen::applyScans(unsigned qId, std::vector<std::shared_ptr<RelationExpr>> &relations) {
    std::shared_ptr<Problem> last = nullptr;
    int id = 0;
    for (auto r : relations) {
        // TODO: Set relation ids in QueryGraph
        r->id = id++;
        std::shared_ptr<Problem> prob = applyScan(qId, r);
        if (last != nullptr)
            last->next = prob;
        else
            dpts[qId].add(prob, 0);
        last = prob;
    }
}

void PlanGen::applyFilters(unsigned qId, std::vector<std::shared_ptr<FilterExpr>> &filters,
                           std::shared_ptr<Problem> prob) {
    for (auto f: filters) { // For each filter in the query
        if (f->applies(prob)) { // Filter applies to problem
            for (auto pl = prob->plans; pl != nullptr; pl = pl->next) { // For each plan
                if (!pl->contains(f)) {
                    auto plan = this->applyFilter(f, pl);
                    this->dpts[qId].add(plan, prob->relations);
                }
            }
        }
    }
}

void PlanGen::applyJoins(unsigned qId, std::vector<std::shared_ptr<JoinExpr>> &joins,
                         std::shared_ptr<Problem> prob1, std::shared_ptr<Problem> prob2) {
    for (auto j: joins) { // For each join
        if (j->applies(prob1, prob2)) { // Apply join and create larger problems
            // For each plan of the 1st problem
            for (auto pl1 = prob1->plans; pl1 != nullptr; pl1 = pl1->next) {
                // For each plan of the 2nd problem
                for (auto pl2 = prob2->plans; pl2 != nullptr; pl2 = pl2->next) {
                    auto plan = j->semi() ? applySemiJoin(j, pl1, pl2) :
                                                            applyJoin(j, pl1, pl2);
                    BitSet m = prob1->relations.unionWith(prob2->relations);
                    dpts[qId].add(plan, m);
                }
            }
        }
    }
}

void PlanGen::applyJoins(unsigned qId, std::vector<std::shared_ptr<JoinExpr>> &joins,
                         std::shared_ptr<Problem> prob) {
    for (auto j: joins) { // For each join
        if (j->applies(prob)) { // Create larger problem
            // For each plan of the problem
            for (auto pl = prob->plans; pl != nullptr; pl = pl->next) {
                if (!pl->contains(j)) {
                    auto plan = applyJoin(j, pl);
                    dpts[qId].add(plan, prob->relations);
                }
            }
        }
    }
}

void PlanGen::applyFunctions(std::shared_ptr<QueryGraph> query, std::shared_ptr<Problem> prob) {
    if (query->tFunctions.empty() || query->groupBy)
        return;
    unsigned min = query->numOp() - query->tFunctions.size();
    unsigned max = query->numOp();
    for (auto pl=prob->plans; pl!=nullptr; pl=pl->next) {
        if ((pl->ops >= min) && (pl->ops < max)) {
            for (auto f : query->tFunctions) {
                if (pl->contains(f))
                    continue;
                auto plan = applyFunction(f, pl);
                dpts[query->id].add(plan, prob->relations);
            }
        }
    }
}

void PlanGen::applyUnaryOperators(unsigned qId, std::shared_ptr<QueryGraph> qGraph,
                                  std::shared_ptr<Problem> prob) {
    // Apply filters
    applyFilters(qId, qGraph->filters, prob);
    // Apply joins that can be reduced into filters
    applyJoins(qId, qGraph->joins, prob);
    // Apply group-by
    applyGroupBy(qGraph, prob);
    // Apply order-by
    applyOrderBy(qGraph, prob);
    // Apply distinct
    applyDistinct(qGraph, prob);
    // Apply table functions in SELECT clause
    applyFunctions(qGraph, prob);
}

std::shared_ptr<Plan> PlanGen::pickOptimalPlan(unsigned qId) {
    std::shared_ptr<Problem> pr = dpts[qId].finalProblem();
    std::shared_ptr<Plan> opt = pr->plans;
    assert(opt != nullptr);
    unsigned s = dpts[qId].stages() - 1;
    for (int i=0; i<dpts[qId].entries(s); i++) {
        auto pr = dpts[qId].at(s, i);
        if (pr == nullptr) continue;
        for (auto pl=pr->plans; pl!=nullptr; pl=pl->next) {
            if (pl->absoluteCost() < opt->absoluteCost()) {
                opt = pl;
            }
        }
    }
    return opt;
}

void PlanGen::applyRules(int qId) {
    // Populate rule base
    std::vector<std::shared_ptr<Rule>> rules;
    rules.push_back(std::make_shared<Rule1>(Rule1()));
    rules.push_back(std::make_shared<Rule2>(Rule2()));
    rules.push_back(std::make_shared<Rule3>(Rule3()));
    rules.push_back(std::make_shared<Rule5>(Rule5()));
    rules.push_back(std::make_shared<Rule6>(Rule6()));
    // Apply rules
    std::shared_ptr<Problem> pr = dpts[qId].finalProblem();
    std::shared_ptr<Plan> p = pr->plans;
    for (; p!=nullptr; p=p->next) {
        for (auto r : rules) {
            if (auto np = r->apply(p)) {
                this->dpts[qId].addPlan(pr, np);
            }
        }
    }
}

// Translates the given query graph into a logical execution plan
std::shared_ptr<Plan> PlanGen::generatePlan(std::shared_ptr<QueryGraph> qGraph,
                                            std::vector<std::shared_ptr<Plan>> &plans,
                                            Database &db) {
    leaderLogging("[INFO] Generating logical plan...\n");
    this->db = &db;
    std::shared_ptr<Plan> opt = nullptr;

    std::map<int, std::shared_ptr<QueryGraph>> qMap;
    qGraph->collectQueries(qMap);
    dpts.resize(qMap.size());
    // Traverse queries and generate plans
    for (auto it=qMap.rbegin(); it!=qMap.rend(); ++it) { // For each sub-query
        int qId = it->first;
        // Initialize plan container
        dpts[qId].allocate(qId, it->second->numOp() + 1, it->second->relations.size());
        // Populate container with scans
        this->applyScans(qId, it->second->relations);
        // Start constructing plans
        for (int i=0; i<dpts[qId].stages(); i++) { // For each stage
            for (int j=0; j<dpts[qId].entries(i); j++) { // For each entry in the current stage
                auto pr1 = dpts[qId].at(i,j);
                for (; pr1!=nullptr; pr1=pr1->next) { // For each problem
                    // Apply unary operators
                    this->applyUnaryOperators(qId, it->second, pr1);
                    // Apply binary operators
                    for (int k=0; k<=i; k++) { // For each previous stage (incl. the current one)
                        for (int m=0; m<dpts[qId].entries(k); m++) { // For each entry of the previous stage
                            auto pr2 = dpts[qId].at(k, m);
                            for (; pr2!=nullptr; pr2=pr2->next) { // For each problem
                                if (pr2->relations.overlapsWith(pr1->relations))
                                    continue;
                                this->applyJoins(qId, it->second->joins, pr1, pr2);
                            }
                        }
                    }
                }
            }
        }
        // Apply transformation rules
        applyRules(qId);

        // Pick optimal plan
        opt = pickOptimalPlan(qId);
        assert(opt);
        plans.push_back(opt);
        applyProject(qId, it->second->projections, opt);
        if (it->second->isSubQuery()) {
            for (auto r: it->second->relAlias) {
                r->cardinality = opt->cardinality;
                r->schema = opt->attributes;
            }
        }
    }
    // Print optimal plans
    if (get_rank()==0) {
        int i = plans.size();
        for (auto p : plans) {
            std::string msg = "\n[INFO] Optimal plan for query " + std::to_string(--i) + " (Cost: " +
                    std::to_string(p->absoluteCost()) + "):\n";
            leaderLogging(msg);
            p->print();
        }
    }
    postProcess(plans);
    return plans.back();
}

void PlanGen::updateAttributes(std::shared_ptr<Plan> op, std::vector<std::shared_ptr<Expression>> &exp) {
    if (exp.empty())
        return;
    for (auto iter=exp.begin(); iter!=exp.end(); ++iter) {
        auto e = *iter;
        if (e->hasType(Expression::DerivedAttributeRef)) {
            auto atts = op->left->attributes;
            if (e->name.find("_BIT") != std::string::npos) {
                assert(op->left);
                for (auto it=atts.rbegin(); it!=atts.rend(); ++it) {
                    auto other = *it;
                    if (other->name.find("[SEL]") != std::string::npos) {
                        *iter = *it;
                        break;
                    }
                }
            }
            if (e->name.find("_AGG") != std::string::npos) {
                for (auto it=atts.rbegin(); it!=atts.rend(); ++it) {
                    auto other = *it;
                    if (other->name.find("SEL") != std::string::npos && other->arithmeticShare) {
                        *iter = *it;
                        break;
                    }
                }
            }
        }
        // Check sub-expressions
        updateAttributes(op, e->expressions);
    }
}

// Refines plan metadata for the code generator
void PlanGen::postProcess(std::vector<std::shared_ptr<Plan>> &plans) {
    // Retrieve share types
    auto rels = db->getRelations();
    std::map<std::string, bool> sTypes;
    for (auto r : rels) {
        for (int i=0; i<r->colNames.size(); i++) {
            bool flag = r->bShares[i] ? true : false;
            sTypes.insert(std::pair<std::string, bool>(r->colNames[i], flag));
        }
    }
    // For all query plans
    for (auto p : plans) {
        // Update plan
        std::vector<std::shared_ptr<Plan>> operators;
        p->collectOperators(operators);
        std::shared_ptr<Plan> prev;
        for (auto o : operators) {
            std::vector<std::shared_ptr<Expression>> a_vec;
            // Update output attributes
            for (auto a : o->attributes) {
                a->setShareTypes(sTypes);
                if (a->hasType(Expression::DerivedAttributeRef)) {
                    auto d = std::dynamic_pointer_cast<DerivedAttributeExpr>(a);
                    assert(d);
                    if (d->name.find("_BIT") != std::string::npos)
                        continue;
                    if (d->name.find("_AGG") != std::string::npos)
                        continue;
                }
                a_vec.push_back(a);
            }
            if (o->op == Plan::NestedLoopJoin) { // Append join attributes
                std::vector<std::shared_ptr<RelationExpr>> leftRel, rightRel;
                assert(o->left && o->right);
                o->left->baseRelations(leftRel);
                o->right->baseRelations(rightRel);
                assert(!leftRel.empty() && !rightRel.empty());
                std::string l, r;
                for (auto rel : leftRel)
                    l += (rel->alias.empty() ? rel->name : rel->alias) + "*";
                l = l.substr(0, l.size()-1);
                for (auto rel : rightRel)
                    r += (rel->alias.empty() ? rel->name : rel->alias) + "*";
                r = r.substr(0, r.size()-1);
                auto j_exp = std::dynamic_pointer_cast<JoinExpr>(o->exp);
                assert(j_exp);
                auto sa = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("SEL",
                                                                                      o->opId));
                sa->arithmeticShare = true;
                auto sb = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("[SEL]",
                                                                                      o->opId));
                sb->arithmeticShare = false;

                if (j_exp->jType == JoinExpr::Semi) {
                    if (!j_exp->agg.empty()) { // Add arithmetic attribute
                        sa->relation = l;
                        a_vec.push_back(sa);
                    }
                }
                else { // Add arithmetic and boolean attributes
                    sa->relation = l + "*" + r;
                    a_vec.push_back(sa);
                    sb->relation = l + "*" + r;
                    a_vec.push_back(sb);
                }
            }
            o->attributes = a_vec;

            // Update logical expressions
            auto exp = o->exp;
            assert(exp);
            switch (exp->eType) {
                case SQLExpression::Scan:
                {
                    auto s = std::dynamic_pointer_cast<ScanExpr>(exp);
                    assert(s);
                    for (auto a : s->rel->schema)
                        a->setShareTypes(sTypes);
                    break;
                }
                case SQLExpression::Filter:
                {
                    auto f = std::dynamic_pointer_cast<FilterExpr>(exp);
                    assert(f);
                    f->selPredicate->setShareTypes(sTypes);
                    for (auto a : f->attributes)
                        a->setShareTypes(sTypes);
                    break;
                }
                case SQLExpression::Join:
                {
                    auto j = std::dynamic_pointer_cast<JoinExpr>(exp);
                    assert(j && j->predicate);
                    j->predicate->setShareTypes(sTypes);
                    for (auto a : j->agg)
                        a->setShareTypes(sTypes);
                    break;
                }
                case SQLExpression::GroupBy:
                {
                    assert(prev);
                    if (prev->op == Plan::NestedLoopJoin) {
                        if (!prev->hasAAttribute()) {
                            auto a = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("SEL",
                                                                                                 prev->opId));
                            a->arithmeticShare = true;
                            a->relation = prev->attributes[0]->relation;
                            prev->attributes.push_back(a);
                        }
                    }
                    auto g = std::dynamic_pointer_cast<GroupByExpr>(exp);
                    assert(g);
                    for (auto a : g->keys)
                        a->setShareTypes(sTypes);
                    updateAttributes(o, g->keys);
                    for (auto a : g->aggs)
                        a->setShareTypes(sTypes);
                    updateAttributes(o, g->aggs);
                    if (g->having)
                        g->having->setShareTypes(sTypes);
                    break;
                }
                case SQLExpression::OrderBy:
                {
                    auto s = std::dynamic_pointer_cast<OrderByExpr>(exp);
                    assert(s);
                    for (auto a : s->keys)
                        a->setShareTypes(sTypes);
                    updateAttributes(o, s->keys);
                    break;
                }
                case SQLExpression::Distinct:
                {
                    auto d = std::dynamic_pointer_cast<DistinctExpr>(exp);
                    assert(d);
                    for (auto a : d->keys)
                        a->setShareTypes(sTypes);
                    updateAttributes(o, d->keys);
                    break;
                }
                case SQLExpression::TableFunction:
                {
                    auto e = std::dynamic_pointer_cast<TableFunctionExpr>(exp);
                    assert(e && e->exp);
                    e->exp->setShareTypes(sTypes);
                    std::vector<std::shared_ptr<Expression>> v = {e->exp};
                    updateAttributes(o, v);
                    break;
                }
                default:
                    leaderLogging("[ERROR] Unrecognized SQL expression\n", true);
            }
            prev = o;
        }
        // Update db metadata
        for (auto r : rels) {
            for (int i=0; i<r->colNames.size(); i++) {
                if (r->bShares[i]) {
                    auto it = r->schema.find(r->colNames[i]);
                    assert(it != r->schema.end());
                    // Update names
                    r->colNames[i] = "[" + r->colNames[i] + "]";
                    // Update index
                    auto c = it->second;
                    r->schema.erase(it);
                    r->schema.insert(std::pair<std::string, Column*>(r->colNames[i], c));
                }
            }
        }
        // leaderLogging("\n[INFO] Transformed query (Cost: :" + std::to_string(p->absoluteCost()) + "):\n");
        // p->print(false, false, false, true);
    }
}