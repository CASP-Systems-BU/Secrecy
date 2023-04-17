#include <iostream>

#include "../../include/planner/plan.h"
#include "../../include/planner/sql_expression.h"


Plan::Plan() {

}

Plan::Plan(const std::shared_ptr<const Plan> other) {
    this->op = other->op;
    this->opId = other->opId;
    this->ordering = other->ordering;
    this->cardinality = other->cardinality;
    this->planCost = other->planCost;
    this->predCost = other->predCost;
    this->attributes = other->attributes;
    this->exp = other->exp;
    this->node = other->node;
    this->nodes = other->nodes;
    this->ops = other->ops;
    this->operators = other->operators;
    this->left = other->left;
    this->right = other->right;
}

Plan::~Plan() {
}

bool Plan::rootIsJoin() const {
    return (op==NestedLoopJoin || op==SemiJoin || op==GroupJoin);
}

bool Plan::hasAAttribute() const {
    for (auto a : attributes) {
        if (a->name.find("SEL") != std::string::npos && a->arithmeticShare)
            return true;
    }
    return false;
}
void Plan::baseRelations(std::vector<std::shared_ptr<RelationExpr>> &relations) const {
    if (op==Scan) {
        auto p = std::dynamic_pointer_cast<ScanExpr>(exp);
        assert(p);
        relations.push_back(p->rel);
    }
    else {
        left->baseRelations(relations);
        if (this->rootIsJoin())
            right->baseRelations(relations);
    }
}

std::shared_ptr<Expression> Plan::getDerivedAttribute() {
    for (auto a : attributes) {
        if (a->hasType(Expression::DerivedAttributeRef))
            return a;
    }
    leaderLogging("Derived attribute not found in plan\n", true);
    return nullptr;
}

void Plan::getKeyPositions(std::vector<unsigned> &key_ops) const {
    if (auto e = std::dynamic_pointer_cast<OrderByExpr>(exp)) {
        for (int i=0; i<e->keys.size(); i++) {
            auto n1 = e->keys[i]->alias.empty() ? e->keys[i]->name : e->keys[i]->alias;
            for (int j=0; j<attributes.size(); j++) {
                auto n2 = attributes[j]->alias.empty() ? attributes[j]->name : attributes[j]->alias;
                if (n1 == n2) {
                    key_ops.push_back(j);
                    break;
                }
            }
        }
    }
}

int Plan::singleOrigin(const std::vector<std::shared_ptr<Expression>> &attributes) const {
    // Fetch left attributes
    std::vector<std::shared_ptr<RelationExpr>> leftRel;
    left->baseRelations(leftRel);
    bool l=false, r=false;
    // Check left input
    for (auto rel : leftRel)
        for (auto a: attributes) {
            if (auto att = std::dynamic_pointer_cast<AttributeExpr>(a))
                l |= rel->contains(att);
            else { // TODO: support arbitrary expressions in DISTINCT
                assert(a->hasType(Expression::DerivedAttributeRef));
            }
        }
    if (this->rootIsJoin()) {
        // Fetch right attributes
        std::vector<std::shared_ptr<RelationExpr>> rightRel;
        right->baseRelations(rightRel);
        // Check right input
        for (auto rel : rightRel)
            for (auto a: attributes) {
                if (auto att = std::dynamic_pointer_cast<AttributeExpr>(a))
                    r |= rel->contains(att);
                else { // TODO: support arbitrary expressions in DISTINCT
                    assert(a->hasType(Expression::DerivedAttributeRef));
                }
            }
    }
    if (l ^ r) {
        if (l) return 0;
        return 1;
    }
    return -1;
}

bool Plan::contains(std::shared_ptr<SQLExpression> op) const {
    assert(op->id != UINT32_MAX);
    return operators.contains(op->id);
}

bool Plan::isSelfJoin(bool adjacent) const {
    if (this->op == JoinFilter)
        return this->left->isSelfJoin(adjacent);
    if ((this->op == NestedLoopJoin) && (this->left->op == Scan) && (this->right->op == Scan)) {
        std::vector<std::shared_ptr<RelationExpr>> leftRel;
        left->baseRelations(leftRel);
        std::vector<std::shared_ptr<RelationExpr>> rightRel;
        right->baseRelations(rightRel);
        assert(!leftRel.empty() && !rightRel.empty());
        if (leftRel[0]->name == rightRel[0]->name) { // Relations are the same
            if (adjacent) { // Check if the join predicate compares adjacent records
                if (auto j = std::dynamic_pointer_cast<JoinExpr>(this->exp)) {
                    auto root = j->findRoot();
                    if (root->isConjunctive()) { // Join predicate is conjunctive
                        // Check if there is a predicate of the form "Equal[r.row_no, Addition[r.row_no, 1]]"
                        std::vector<std::shared_ptr<Expression>> atoms;
                        root->flatten(atoms);
                        for (auto a : atoms) {
                            if (a->hasType(Expression::Equal)) {
                                auto l = a->expressions[0];
                                auto r = a->expressions[1];
                                if (l->name=="row_no") {
                                    assert(r->hasType(Expression::Addition));
                                    if (r->expressions[0]->name=="row_no" && r->expressions[1]->ival==1) {
                                        a->adjRows = true;
                                        return true;
                                    }
                                }
                            }
                        }
                    }
                    return false;
                }
                leaderLogging("[ERROR] Unrecognized self-join expression\n", true);
            }
            return true;
        }
    }
    return false;
}

bool Plan::isJoin(std::vector<std::shared_ptr<SQLExpression>> &predicates) const {
    if (this->op == JoinFilter) {
        predicates.push_back(this->exp);
        return this->left->isJoin(predicates);
    }
    if (this->op == NestedLoopJoin) {
        predicates.push_back(this->exp);
        std::vector<std::shared_ptr<RelationExpr>> leftRel;
        left->baseRelations(leftRel);
        std::vector<std::shared_ptr<RelationExpr>> rightRel;
        right->baseRelations(rightRel);
        assert(!leftRel.empty() && !rightRel.empty());
        if (leftRel.size() != rightRel.size())
            return true;
        for (int i; i<leftRel.size(); i++) {
            if (leftRel[i]->name != rightRel[i]->name)
                return true;
        }
    }
    return false;
}

void Plan::collectScans(std::vector<std::shared_ptr<Plan>> &scans) {
    if (op == Scan)
        scans.push_back(shared_from_this());
    else {
        assert(left != nullptr);
        left->collectScans(scans);
        if (right)
            right->collectScans(scans);
    }
}

void Plan::collectOps(std::vector<std::shared_ptr<Plan>> &operators) {
    if (left)
        operators.push_back(left);
    if (right)
        operators.push_back(right);
    if (left)
        left->collectOps(operators);
    if (right)
        right->collectOps(operators);
}

void Plan::collectOperators(std::vector<std::shared_ptr<Plan>> &operators) {
    operators.push_back(shared_from_this());
    collectOps(operators);
    std::reverse(operators.begin(), operators.end());
}

void Plan::getJoinInputs(std::vector<std::shared_ptr<Plan>> &plans) {
    if (op == JoinFilter) {
        this->left->getJoinInputs(plans);
    }
    else if (op == NestedLoopJoin) {
        plans.push_back(left);
        plans.push_back(right);
    }
}

bool Plan::isDistinct() const {
    assert(exp != nullptr);
    if (auto d = std::dynamic_pointer_cast<DistinctExpr>(exp))
        return true;
    return false;
}

// Checks if `this` plan dominates another plan
bool Plan::dominates(const std::shared_ptr<Plan> other) const {
    if (this->ordering == other->ordering &&
        this->operators == other->operators &&
        this->planCost < other->planCost)
        return true;
    return false;
}

// Returns the absolute plan cost based on parameters ALPHA and BETA
unsigned Plan::absoluteCost() const {
   return ALPHA * (this->planCost.operation + this->planCost.compOp)
        + BETA * (this->planCost.synchronization + this->planCost.compSync);
}

// Prints the plan to stdout
void Plan::print(unsigned indent, bool pCost, bool rootOnly, bool schema) const {
    if(get_rank() == 0) {
        assert(exp != nullptr);
        std::string ind;
        for (unsigned index = 0; index < indent; index++)
            ind += ' ';
        if (pCost) {
            std::cout << ind << "Cost: " << absoluteCost() << std::endl;
        }
        if (schema) {
            std::string sch = ind + "[Schema: ";
            for (auto a : attributes)
                sch += a->toString() + ", ";
            if (attributes.size() > 0)
                sch = sch.substr(0, sch.size() - 2);
            sch += "]";
            std::cout << sch << std::endl;
        }
        std::cout << ind;
        switch (op) {
            case Scan: {
                auto s = std::dynamic_pointer_cast<ScanExpr>(this->exp);
                std::cout << "Scan (" << s->rel->toString() << ")" << std::endl;
                break;
            }
            case NestedLoopJoin:
                std::cout << "NestedLoopJoin [Type: " << exp->toString() << "]" << std::endl;
                break;
            case Filter:
            case SemiJoin:
            case GroupBy:
            case OrderBy:
            case Distinct:
            case TableFunction:
            case AdjacentEq:
            case OddEvenMerge:
            case Mask:
                std::cout << exp->toString() << std::endl;
                break;
            case JoinFilter:
                assert(this->exp);
                std::cout << "JoinFilter [Type: " << this->exp->toString() << "]" << std::endl;
                break;
            case Project:
            case AdjacentGeq:
            case Composition:
            case RowFunction:
            case Union:
            case GroupJoin:
                break;
        }
        if (!rootOnly) {
            // Continue recursively
            switch (op) {
                case Scan:
                    break;
                case NestedLoopJoin:
                case SemiJoin:
                case Union:
                    assert(left != nullptr && right != nullptr);
                    left->print(indent + 1, pCost, rootOnly, schema);
                    right->print(indent + 1, pCost, rootOnly, schema);
                    break;
                case Project:
                case Filter:
                case JoinFilter:
                case TableFunction:
                case GroupBy:
                case OrderBy:
                case Distinct:
                case AdjacentEq:
                case OddEvenMerge:
                case Mask:
                    assert(left != nullptr);
                    left->print(indent + 1, pCost, rootOnly, schema);
                    break;
                case GroupJoin:
                case AdjacentGeq:
                case Composition:
                case RowFunction:
                    break;
            }
        }
    }
}