#include "../../include/planner/rule.h"

std::shared_ptr<Plan> Rule::apply(std::shared_ptr<Plan> plan) {
    if (auto info = match(plan))
        return transform(plan, info);
    return nullptr;
}

std::shared_ptr<Rule::MatchInfo> Rule1::match(std::shared_ptr<Plan> plan) {
    auto distinct = plan;
    if (distinct->op == Plan::Op::Distinct) {
        auto join = distinct->left;
        assert(join != nullptr);
        if (join->op == Plan::NestedLoopJoin) {
            auto origin = join->singleOrigin(distinct->attributes);
            if (origin >= 0) { // We have a match
                auto info = std::make_shared<Rule::MatchInfo>(Rule::MatchInfo());
                info->op1 = distinct;
                info->op2 = join;
                info->origin = origin;
                return info;
            }
        }
    }
    return nullptr;
}

std::shared_ptr<Plan> Rule1::transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info) {
    auto distinct = info->op1;
    auto join = info->op2;
    assert(distinct && join);
    std::shared_ptr<Plan> l, r;
    if (info->origin == 0) {
        l = join->left;
        r = join->right;
    }
    else {
        assert(info->origin == 1);
        l = join->right;
        r = join->left;
    }
    // Step 0: Create semi-join
    auto semi = std::make_shared<Plan>(Plan(join));
    semi->op = Plan::NestedLoopJoin;
    // Create semi-join expression
    auto jExpr = std::dynamic_pointer_cast<JoinExpr>(join->exp);
    auto sjExpr = std::make_shared<JoinExpr>(JoinExpr(JoinExpr::Semi,
                                                      jExpr->predicate, jExpr->left,
                                                      jExpr->right, join->opId));
    semi->exp = sjExpr;
    semi->attributes = l->attributes;
    auto e = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("SJ_BIT",
                                                                                        sjExpr->id));
    semi->attributes.push_back(e);
    semi->left = l;
    semi->right = r;
    semi->cardinality = l->cardinality;
    semi->planCost = CostModel::semiJoinCost(l->cardinality, r->cardinality,
                                             semi->predCost);
    // Step 1: Add sort on top of semi-join
    auto sort = std::make_shared<Plan>(Plan());
    sort->op = Plan::OrderBy;
    sort->opId = distinct->opId;
    auto oExpr = std::make_shared<OrderByExpr>(OrderByExpr(distinct->opId));
    auto d = std::dynamic_pointer_cast<DistinctExpr>(distinct->exp);
    assert(d);
    oExpr->keys.resize(d->keys.size() + 1);
    oExpr->asc.resize(d->keys.size() + 1, true);
    oExpr->keys[0] = e; // First sorting key
    for (int i = 0; i < d->keys.size(); i++)
        oExpr->keys[i + 1] = d->keys[i];  // Remaining sorting keys
    sort->exp = oExpr;
    sort->nodes = semi->nodes;
    sort->operators = semi->operators;
    sort->ops = semi->ops + 1;
    sort->cardinality = semi->cardinality;
    sort->left = semi;
    sort->attributes = semi->attributes;
    sort->planCost = sort->left->planCost + CostModel::sortCost(sort->cardinality);
    // Step 2: Add adjacent eq on top of sort
    auto eq = std::make_shared<Plan>(Plan());
    eq->op = Plan::AdjacentEq;
    eq->opId = distinct->opId;
    auto eqExpr = std::make_shared<DistinctExpr>(DistinctExpr(distinct->opId));
    eqExpr->adjacentEq = true;
    eqExpr->keys = d->keys;
    eq->exp = eqExpr;
    for (auto a : sort->attributes)
        if (!a->hasType(Expression::DerivedAttributeRef))
            eq->attributes.push_back(a);
    eq->nodes = sort->nodes;
    eq->ops = sort->ops;
    eq->operators = sort->operators;
    eq->cardinality = sort->cardinality;
    eq->left = sort;
    auto dbit = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("D_BIT"));
    eq->attributes.push_back(dbit);
    eq->planCost = eq->left->planCost + CostModel::adjacentEq(eq->cardinality);
    // Step 3: Add mask
    auto m = std::make_shared<Plan>(Plan());
    m->op = Plan::Mask;
    m->opId = eq->opId;
    auto expr = std::make_shared<Expression>(Expression(Expression::Equal));
    auto c = std::make_shared<ConstantExpr>(ConstantExpr((int64_t) 0));
    expr->expressions.push_back(dbit);
    expr->expressions.push_back(c);
    auto mExpr = std::make_shared<TableFunctionExpr>(TableFunctionExpr(expr, distinct->opId));
    mExpr->mask = true;
    m->exp = mExpr;
    m->nodes = eq->nodes;
    m->ops = eq->ops;
    m->operators = eq->operators;
    m->cardinality = eq->cardinality;
    m->attributes = eq->attributes;
    m->left = eq;
    auto eqC = CostModel::eqCost();
    m->planCost = m->left->planCost + CostModel::selectCost(m->cardinality, eqC);
    return m;
}

std::shared_ptr<Rule::MatchInfo> Rule2::match(std::shared_ptr<Plan> plan) {
    auto group = plan;
    if (group->op == Plan::Op::GroupBy) {
        auto g = std::dynamic_pointer_cast<GroupByExpr>(group->exp);
        if (g->aggs.size() == 1 && g->aggs[0]->name == "COUNT") {
            auto join = group->left;
            assert(join != nullptr);
            if (join->op == Plan::NestedLoopJoin) {
                auto g = std::dynamic_pointer_cast<GroupByExpr>(group->exp);
                auto origin = join->singleOrigin(g->keys);
                if (origin >= 0) { // We have a match
                    auto info = std::make_shared<Rule::MatchInfo>(Rule::MatchInfo());
                    info->op1 = group;
                    info->op2 = join;
                    info->origin = origin;
                    return info;
                }
            }
        }
    }
    else if (group->left){
        return this->match(group->left);
    }
    return nullptr;
}

std::shared_ptr<Plan> Rule2::transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info) {
    auto group = info->op1;
    auto join = info->op2;
    assert(group && join);
    auto g = std::dynamic_pointer_cast<GroupByExpr>(group->exp);
    assert(g);
    // Step 0: create semi-join with partial aggregation
    auto semi = std::make_shared<Plan>(Plan(join));
    auto j = std::dynamic_pointer_cast<JoinExpr>(join->exp);
    std::shared_ptr<FunctionExpr> cnt;
    if (j->jType == JoinExpr::Inner || j->jType == JoinExpr::LeftOuter) { // Add partial aggregation
        cnt = std::make_shared<FunctionExpr>(FunctionExpr("COUNT", nullptr));
        auto att = std::make_shared<Expression>(Expression(Expression::Star));
        cnt->expressions.push_back(att);
    }
    semi->op = Plan::NestedLoopJoin;
    auto jExpr = std::dynamic_pointer_cast<JoinExpr>(join->exp);
    auto sjExpr = std::make_shared<JoinExpr>(JoinExpr(JoinExpr::Semi,
                                                      jExpr->predicate, jExpr->left,
                                                      jExpr->right, join->opId));
    if (j->jType == JoinExpr::Inner || j->jType == JoinExpr::LeftOuter) {
        sjExpr->agg.push_back(cnt);
    }
    std::shared_ptr<Plan> l, r;
    if (info->origin == 0) {
        l = join->left;
        r = join->right;
    }
    else {
        assert(info->origin == 1);
        l = join->right;
        r = join->left;
    }
    semi->exp = sjExpr;
    semi->planCost = join->left->planCost + CostModel::semiJoinCost(l->cardinality,
                                                                    r->cardinality,
                                                                    semi->predCost);
    semi->attributes = l->attributes;
    auto a = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("SJ_AGG"));
    if (!sjExpr->agg.empty()) {
        semi->attributes.push_back(a);
    }
    // Step 1: Sort left semi-join input
    auto sort = std::make_shared<Plan>(Plan());
    sort->op = Plan::OrderBy;
    sort->opId = group->opId;
    auto oExpr = std::make_shared<OrderByExpr>(OrderByExpr(group->opId));
    oExpr->keys = g->keys;
    oExpr->asc.resize(g->keys.size(), true);
    sort->exp = oExpr;
    sort->left = l;
    sort->nodes = l->nodes;
    sort->operators = l->operators;
    sort->ops = l->ops + 1;
    sort->cardinality = l->cardinality;
    sort->attributes = l->attributes;
    sort->planCost = l->planCost + CostModel::sortCost(sort->cardinality);
    semi->left = sort;
    semi->right = r;
    // Step 2: Add odd-even aggregation
    auto oddEven = std::make_shared<Plan>(Plan());
    oddEven->op = Plan::OddEvenMerge;
    oddEven->opId = group->opId;
    auto gAgg = std::make_shared<FunctionExpr>(FunctionExpr("SUM", nullptr));
    gAgg->expressions.push_back(a);
    auto oeExpr = std::make_shared<GroupByExpr>(GroupByExpr(group->opId));
    oeExpr->oddEven = true;
    oeExpr->keys = g->keys;
    oeExpr->aggs.push_back(gAgg);
    oddEven->exp = oeExpr;
    oddEven->left = semi;
    oddEven->nodes = semi->nodes;
    oddEven->operators = group->operators;
    oddEven->ops = semi->ops + 1;
    oddEven->cardinality = semi->cardinality;
    oddEven->attributes = semi->attributes;
    oddEven->attributes.push_back(gAgg);
    // Aggregation function is local after single-bit conversion
    auto c = CostModel::Cost(1,0);
    oddEven->planCost = oddEven->left->planCost +
                        CostModel::oddEvenMerge(oddEven->cardinality, c) +
                        CostModel::bitConvCost(oddEven->cardinality);
    return oddEven;
}

std::shared_ptr<Rule::MatchInfo> Rule3::match(std::shared_ptr<Plan> plan) {
    auto distinct = plan;
    if (distinct->op == Plan::Op::Distinct) {
        auto select = distinct->left;
        assert(select != nullptr);
        if (select->op == Plan::Filter) { // We have a match
            auto info = std::make_shared<Rule::MatchInfo>(Rule::MatchInfo());
            info->op1 = distinct;
            info->op2 = select;
            return info;
        }
    }
    return nullptr;
}

std::shared_ptr<Plan> Rule3::transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info) {
    auto distinct = info->op1;
    auto select = info->op2;
    // Step 0: Sort on a_phi, distinct keys
    auto sort = std::make_shared<Plan>(Plan());
    sort->op = Plan::OrderBy;
    sort->opId = distinct->opId;
    // Create order-by expression
    auto d = std::dynamic_pointer_cast<DistinctExpr>(distinct->exp);
    assert(d);
    auto oExpr = std::make_shared<OrderByExpr>(OrderByExpr(distinct->opId));
    oExpr->asc.resize(d->keys.size() + 1, true);
    auto a = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("S_BIT"));
    oExpr->keys.push_back(a);
    oExpr->keys.insert(oExpr->keys.end(), d->keys.begin(), d->keys.end());
    sort->exp = oExpr;
    sort->left = select;
    sort->nodes = select->nodes;
    sort->operators = select->operators;
    sort->ops = select->ops + 1;
    sort->cardinality = select->cardinality;
    sort->attributes = select->attributes;
    sort->planCost = sort->left->planCost + CostModel::sortCost(sort->cardinality);
    // Step 1: Add adjacent geq on distinct keys
    auto eq = std::make_shared<Plan>(Plan());
    eq->op = Plan::AdjacentEq;
    eq->opId = distinct->opId;
    auto eqExpr = std::make_shared<DistinctExpr>(DistinctExpr(distinct->opId));
    eqExpr->adjacentEq = true;
    eqExpr->keys = d->keys;
    eq->exp = eqExpr;
    eq->nodes = sort->nodes;
    eq->ops = sort->ops;
    eq->operators = sort->operators;
    eq->cardinality = sort->cardinality;
    eq->left = sort;
    auto dbit = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("D_BIT"));
    for (auto a : sort->attributes)
        if (!a->hasType(Expression::DerivedAttributeRef))
            eq->attributes.push_back(a);
    eq->attributes.push_back(dbit);
    eq->planCost = eq->left->planCost + CostModel::adjacentEq(eq->cardinality);
    // Step 2: Add mask
    auto m = std::make_shared<Plan>(Plan());
    m->op = Plan::Mask;
    m->opId = eq->opId;
    m->ops = eq->ops;
    auto expr = std::make_shared<Expression>(Expression(Expression::Equal));
    auto c = std::make_shared<ConstantExpr>(ConstantExpr((int64_t) 0));
    expr->expressions.push_back(dbit);
    expr->expressions.push_back(c);
    auto mExpr = std::make_shared<TableFunctionExpr>(TableFunctionExpr(expr, distinct->opId));
    mExpr->mask = true;
    m->exp = mExpr;
    m->nodes = eq->nodes;
    m->operators = eq->operators;
    m->cardinality = eq->cardinality;
    m->left = eq;
    auto eqC = CostModel::eqCost();
    m->attributes = eq->attributes;
    m->planCost = m->left->planCost + CostModel::selectCost(m->cardinality, eqC);
    return m;
}

std::shared_ptr<Rule::MatchInfo> Rule5::match(std::shared_ptr<Plan> plan) {
    auto distinct = plan;
    if (distinct->isDistinct()) {
        auto join = distinct->left;
        assert(join != nullptr);
        if (join->isSelfJoin(true)) { // We have a match
            auto info = std::make_shared<Rule::MatchInfo>(Rule::MatchInfo());
            info->op1 = distinct;
            info->op2 = join;
            return info;
        }
    }
    return nullptr;
}

std::shared_ptr<Plan> Rule5::transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info) {
    auto distinct = info->op1;
    auto join = info->op2;
    assert(distinct && join);
    auto j = std::dynamic_pointer_cast<JoinExpr>(join->exp);
    assert(j);
    auto root = j->predicate->findRoot();
    assert(root);
    std::vector<std::shared_ptr<Expression>> atoms, eq_atoms, rem_atoms;
    root->flatten(atoms);
    for (auto a : atoms) {
        if (a->hasType(Expression::Equal)) {
            auto l = a->expressions[0];
            auto r = a->expressions[1];
            if (a->adjRows)
                continue;
            else if (l->hasType(Expression::AttributeRef) &&
                     r->hasType(Expression::AttributeRef)) {
                eq_atoms.push_back(a);
                continue;
            }
        }
        rem_atoms.push_back(a);
    }
    // Step 0: Sort on attributes in equality predicates
    auto sort = std::make_shared<Plan>(Plan());
    sort->op = Plan::OrderBy;
    sort->opId = join->opId;
    // Create order-by expression
    std::vector<std::shared_ptr<Plan>> scans;
    join->collectScans(scans);
    assert(!scans.empty());
    auto s = std::dynamic_pointer_cast<ScanExpr>(scans[0]->exp);
    assert(s);
    auto oExpr = std::make_shared<OrderByExpr>(OrderByExpr(join->opId));
    oExpr->asc.resize(eq_atoms.size(), true);
    for (auto a : eq_atoms) {
        if (s->rel->alias == a->expressions[0]->relation)
            oExpr->keys.push_back(a->expressions[0]);
        else if (s->rel->alias == a->expressions[1]->relation)
            oExpr->keys.push_back(a->expressions[1]);
    }
    sort->exp = oExpr;
    sort->left = scans[0];
    sort->nodes = 1;
    sort->operators.set(sort->opId);
    sort->ops = 1;
    sort->cardinality = scans[0]->cardinality;
    sort->attributes = scans[0]->attributes;
    sort->planCost = CostModel::sortCost(sort->cardinality);
    // Step 1: Apply join predicates
    std::vector<std::shared_ptr<SQLExpression>> predicates;
    for (auto a : eq_atoms) {
        auto p = std::make_shared<TableFunctionExpr>(TableFunctionExpr(a, a->opId));
        p->adjacent = true;
        predicates.push_back(p);
    }
    for (auto a : rem_atoms) {
        auto p = std::make_shared<TableFunctionExpr>(TableFunctionExpr(a, a->opId));
        p->adjacent = true;
        predicates.push_back(p);
    }
    std::shared_ptr<Plan> prev = nullptr;
    for (auto p : predicates) {
        auto tf = std::make_shared<Plan>(Plan());
        tf->op = Plan::TableFunction;
        tf->opId = p->id;
        tf->exp = p;
        if (prev == nullptr) {
            tf->left = sort;
            tf->attributes = sort->attributes;
        }
        else {
            tf->left = prev;
            for (auto a : prev->attributes)
                if (!a->hasType(Expression::DerivedAttributeRef))
                    tf->attributes.push_back(a);
        }
        tf->attributes.push_back(std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("J_BIT",
                                                                                             tf->opId)));
        prev = tf;
        tf->nodes = tf->left->nodes + 1;
        tf->operators.set(tf->opId);
        tf->ops = tf->left->ops + 1;
        tf->cardinality = tf->left->cardinality;
        tf->planCost = CostModel::ineqCost(tf->cardinality);

    }
    // Step 2: Sort on join_bit, distinct key
    auto sort2 = std::make_shared<Plan>(Plan());
    sort2->op = Plan::OrderBy;
    sort2->opId = distinct->opId;
    // Create order-by expression
    auto oExpr2 = std::make_shared<OrderByExpr>(OrderByExpr(distinct->opId));
    auto d = std::dynamic_pointer_cast<DistinctExpr>(distinct->exp);
    oExpr2->asc.resize(d->keys.size() + 1, true);
    auto jBit = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("J_BIT",
                                                                                            prev->opId));
    oExpr2->keys.push_back(jBit);
    oExpr2->keys.insert(oExpr2->keys.end(), d->keys.begin(), d->keys.end());
    sort2->exp = oExpr2;
    sort2->left = prev;
    sort2->nodes = prev->nodes;
    sort2->operators.set(sort2->opId);
    sort2->ops = prev->ops + 1;
    sort2->cardinality = prev->cardinality;
    sort2->planCost = CostModel::sortCost(sort2->cardinality);
    sort2->attributes = prev->attributes;
    // Step 3: Add adjacent eq on distinct key
    auto eq = std::make_shared<Plan>(Plan());
    eq->op = Plan::AdjacentEq;
    eq->opId = distinct->opId;
    for (auto a : prev->attributes)
        if (!a->hasType(Expression::DerivedAttributeRef))
            eq->attributes.push_back(a);
    auto dbit = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("D_BIT"));
    eq->attributes.push_back(dbit);
    auto eqExpr = std::make_shared<DistinctExpr>(DistinctExpr(distinct->opId));
    eqExpr->adjacentEq = true;
    eqExpr->keys = d->keys;
    eq->exp = eqExpr;
    eq->nodes = sort2->nodes;
    eq->ops = sort2->ops;
    eq->operators = sort2->operators;
    eq->cardinality = sort2->cardinality;
    eq->left = sort2;
    eq->planCost = eq->left->planCost + CostModel::adjacentEq(eq->cardinality);
    // Step 4: Add masking
    auto m = std::make_shared<Plan>(Plan());
    m->op = Plan::Mask;
    m->opId = eq->opId;
    m->ops = eq->ops;
    auto expr = std::make_shared<Expression>(Expression(Expression::Equal));
    auto c = std::make_shared<ConstantExpr>(ConstantExpr((int64_t) 0));
    expr->expressions.push_back(dbit);
    expr->expressions.push_back(c);
    auto mExpr = std::make_shared<TableFunctionExpr>(TableFunctionExpr(expr,
                                                                                    distinct->opId));
    mExpr->mask = true;
    m->exp = mExpr;
    m->nodes = eq->nodes;
    m->operators = eq->operators;
    m->cardinality = eq->cardinality;
    m->left = eq;
    auto eqC = CostModel::eqCost();
    m->planCost = m->left->planCost + CostModel::selectCost(m->cardinality, eqC);
    m->attributes = eq->attributes;
    return m;
}

std::shared_ptr<Rule::MatchInfo> Rule6::match(std::shared_ptr<Plan> plan) {
    auto tf = plan;
    if (auto f = std::dynamic_pointer_cast<TableFunctionExpr>(tf->exp)) {
        if (f->exp->name == "COUNT" && f->exp->distinct) {
            auto join = tf->left;
            assert(join != nullptr);
            std::vector<std::shared_ptr<SQLExpression>> predicates;
            if (join->isJoin(predicates)) {
                // Check if join predicate has the form And[Equal[r.dist_key s.dist_key], LessThanEqual[r.a, s.b]]
                auto j = std::dynamic_pointer_cast<JoinExpr>(join->exp);
                auto root = j->findRoot();
                if (root->isConjunctive() && predicates.size() == 2) {
                    std::shared_ptr<Expression> l = nullptr, r = nullptr;
                    auto p1 = std::dynamic_pointer_cast<JoinExpr>(predicates[0]);
                    auto p2 = std::dynamic_pointer_cast<JoinExpr>(predicates[1]);
                    assert(p1 && p2);
                    if (p1->predicate->hasType(Expression::Equal)) {
                        l = p1->predicate;
                        r = p2->predicate;
                    } else {
                        l = p2->predicate;
                        r = p1->predicate;
                    }
                    if (l->hasType(Expression::Equal)) { // 1st predicate is equality
                        auto ll = l->expressions[0];
                        auto lr = l->expressions[1];
                        if (ll->hasType(Expression::AttributeRef) &&
                            lr->hasType(Expression::AttributeRef)) { // Both are attributes
                            auto k = f->exp->expressions[0];
                            if (k->hasType(Expression::AttributeRef) &&
                                k->name == ll->name &&
                                k->name == lr->name) { // Both attributes in the equality have the same name as the
                                // distinct attribute
                                if (r->hasType(Expression::LessThanEqual)) {// 2nd predicate is inequality
                                    auto rl = r->expressions[0];
                                    auto rr = r->expressions[1];
                                    if (rl->hasType(Expression::AttributeRef) &&
                                        rr->hasType(Expression::AttributeRef)) { // Both are attributes
                                        if (rl->relation == ll->relation &&
                                            rr->relation == lr->relation) { // We have a match
                                            auto info = std::make_shared<MatchInfo>(MatchInfo());
                                            info->op1 = tf;
                                            info->op2 = join;
                                            info->predicates = predicates;
                                            info->e1 = rl;
                                            info->e2 = rr;
                                            return info;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return nullptr;
}

std::shared_ptr<Plan> Rule6::transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info) {
    auto tf = info->op1;
    auto join = info->op2;
    auto rl = info->e1;
    auto rr = info->e2;
    assert(tf && join && rl && rr);
    auto f = std::dynamic_pointer_cast<TableFunctionExpr>(tf->exp);
    assert(f);
    std::vector<std::shared_ptr<Plan>> plans;
    join->getJoinInputs(plans);
    assert(plans.size()==2);
    std::vector<std::shared_ptr<RelationExpr>> lRels, rRels;
    plans[0]->baseRelations(lRels);
    plans[1]->baseRelations(rRels);
    assert(lRels.size()==1 && rRels.size()==1);
    // Step 0: Apply left sort
    auto sort = std::make_shared<Plan>(Plan());
    sort->op = Plan::OrderBy;
    sort->opId = tf->opId;
    sort->attributes = plans[0]->attributes;
    unsigned ks = plans[0]->requiresComposition() ? f->exp->expressions.size() + 2 :
                                                    f->exp->expressions.size() +1;
    // Create order-by expression
    auto oExpr = std::make_shared<OrderByExpr>(OrderByExpr(tf->opId));
    oExpr->asc.resize(ks, true);
    if (plans[0]->requiresComposition())
        oExpr->keys.push_back(plans[0]->getDerivedAttribute());
    auto key = std::dynamic_pointer_cast<AttributeExpr>(f->exp->expressions[0]);
    assert(key);
    // Use the right keys and sort direction
    for (auto a : lRels[0]->schema)
        if (key->name == a->name) {
            oExpr->keys.push_back(a);
        }
    if (lRels[0]->identifiedBy(rl->relation)) {
        oExpr->keys.push_back(rl);
    }
    else {
        assert(lRels[0]->identifiedBy(rr->relation));
        oExpr->keys.push_back(rr);
        oExpr->asc[oExpr->asc.size()-1] = false;
    }
    sort->exp = oExpr;
    sort->left = plans[0];
    sort->nodes = sort->left->nodes;
    sort->operators.set(sort->opId);
    sort->ops = plans[0]->nodes + 1;
    sort->cardinality = plans[0]->cardinality;
    sort->planCost = CostModel::sortCost(sort->cardinality);
    // Step 1: Apply right sort
    auto sort2 = std::make_shared<Plan>(Plan());
    sort2->op = Plan::OrderBy;
    sort2->opId = tf->opId;
    sort2->attributes = plans[1]->attributes;
    unsigned ks2 = plans[1]->requiresComposition() ? f->exp->expressions.size() + 2 :
                                                     f->exp->expressions.size() + 1;
    // Create order-by expression
    auto oExpr2 = std::make_shared<OrderByExpr>(OrderByExpr(tf->opId));
    oExpr2->asc.resize(ks2, true);
    if (plans[1]->requiresComposition()) {
        oExpr2->keys.push_back(plans[0]->getDerivedAttribute());
    }
    // Use the right keys and sort direction
    for (auto a : rRels[0]->schema)
        if (key->name == a->name) {
            oExpr2->keys.push_back(a);
        }
    if (rRels[0]->identifiedBy(rr->relation)) {
        oExpr2->keys.push_back(rr);
        oExpr2->asc[oExpr2->asc.size()-1] = false;
    }
    else {
        assert(rRels[0]->identifiedBy(rl->relation));
        oExpr2->keys.push_back(rl);
    }
    sort2->exp = oExpr2;
    sort2->left = plans[1];
    sort2->nodes = sort2->left->nodes;
    sort2->operators.set(sort2->opId);
    sort2->ops = plans[1]->nodes + 1;
    sort2->cardinality = plans[1]->cardinality;
    sort2->planCost = CostModel::sortCost(sort2->cardinality);
    // Step 2: Apply left distinct
    auto eq = std::make_shared<Plan>(Plan());
    eq->op = Plan::AdjacentEq;
    eq->opId = tf->opId;
    for (auto a: sort->attributes)
        if (!a->hasType(Expression::DerivedAttributeRef))
            eq->attributes.push_back(a);
    eq->attributes.push_back(std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("D_BIT",
                                                                                         eq->opId)));
    auto eqExpr = std::make_shared<DistinctExpr>(DistinctExpr(tf->opId));
    eqExpr->adjacentEq = true;
    eqExpr->keys = f->exp->expressions;
    eq->exp = eqExpr;
    eq->nodes = sort->nodes;
    eq->ops = sort->ops;
    eq->operators = sort->operators;
    eq->cardinality = sort->cardinality;
    eq->left = sort;
    eq->planCost = eq->left->planCost + CostModel::adjacentEq(eq->cardinality);
    // Step 3: Apply right distinct
    auto eq2 = std::make_shared<Plan>(Plan());
    eq2->op = Plan::AdjacentEq;
    eq2->opId = tf->opId;
    for (auto a: sort2->attributes)
        if (!a->hasType(Expression::DerivedAttributeRef))
            eq2->attributes.push_back(a);
    eq2->attributes.push_back(std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("D_BIT",
                                                                                          eq2->opId)));
    auto eqExpr2 = std::make_shared<DistinctExpr>(DistinctExpr(tf->opId));
    eqExpr2->adjacentEq = true;
    eqExpr2->keys = f->exp->expressions;
    eq2->exp = eqExpr2;
    eq2->nodes = sort2->nodes;
    eq2->ops = sort2->ops;
    eq2->operators = sort2->operators;
    eq2->cardinality = sort2->cardinality;
    eq2->left = sort2;
    eq2->planCost = eq2->left->planCost + CostModel::adjacentEq(eq2->cardinality);
    // Step 2: Apply join
    std::shared_ptr<Plan> prev;
    for (auto p : info->predicates) {
        auto join2 = std::make_shared<Plan>(Plan());
        join2->opId = join->opId;
        auto joinExpr = std::dynamic_pointer_cast<JoinExpr>(p);
        assert(joinExpr);
        join2->exp = joinExpr;
        if (!prev) {
            join2->op = Plan::NestedLoopJoin;
            join2->nodes = eq->nodes + eq2->nodes;
            join2->ops = eq->ops + eq2->ops;
            join2->operators = eq->operators.unionWith(eq2->operators);
            join2->cardinality = eq->cardinality * eq2->cardinality;
            join2->left = eq;
            join2->right = eq2;
            CostModel::Cost eqC = CostModel::eqCost();
            join2->planCost = eq->planCost + eq2->planCost +
                              CostModel::joinCost(eq->cardinality,
                                                  eq2->cardinality, eqC);
            for (auto a : eq->attributes)
                if (!a->hasType(Expression::DerivedAttributeRef))
                    join2->attributes.push_back(a);
            for (auto a : eq2->attributes)
                if (!a->hasType(Expression::DerivedAttributeRef))
                    join2->attributes.push_back(a);
        }
        else {
            join2->op = Plan::JoinFilter;
            join2->nodes = prev->nodes;
            join2->ops = prev->ops + 1;
            join2->operators = prev->operators;
            join2->cardinality = prev->cardinality;
            join2->left = prev;
            join2->planCost = prev->planCost +
                              CostModel::eqCost(prev->cardinality);
            for (auto a : prev->attributes)
                if (!a->hasType(Expression::DerivedAttributeRef))
                    join2->attributes.push_back(a);
        }
        join2->attributes.push_back(std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("J_BIT",
                                                                                                join2->opId)));
        prev = join2;
    }
    // Step 3: Apply table function
    auto fExp = std::make_shared<TableFunctionExpr>(TableFunctionExpr(f));
    auto a = std::make_shared<DerivedAttributeExpr>(DerivedAttributeExpr("J_BIT"));
    fExp->exp->distinct = false;
    fExp->exp->expressions.clear();
    fExp->exp->expressions.push_back(a);
    auto tf2 = std::make_shared<Plan>(Plan());
    tf2->op = Plan::TableFunction;
    tf2->opId = tf->opId;
    tf2->exp = fExp;
    tf2->nodes = prev->nodes;
    prev->ops = prev->ops + 1;
    prev->operators = prev->operators;
    tf2->cardinality = tf2->cardinality;
    tf2->attributes.push_back(fExp->exp);
    tf2->left = prev;
    return tf2;
}
