#include "../../include/planner/queryGraph.h"
#include <set>

//using namespace std;


void QueryGraph::getRelations(const std::string &rel_id, std::vector<std::shared_ptr<RelationExpr> > &v,
                                     bool where) {
    for (auto r : this->relations) {
        if ((r->name == rel_id) ||
            (!r->alias.empty() && (r->alias == rel_id)))
            v.push_back(r);
    }
    if (where && (whereClause != nullptr)) { // Check if relation exists in WHERE clause
        auto r = whereClause->getRelation(rel_id);
        if (r != nullptr) {
            v.push_back(std::dynamic_pointer_cast<RelationExpr>(r));
        }
    }
}

std::shared_ptr<RelationExpr> QueryGraph::getRelation(const std::string &rel_id, bool where) const {
    for (auto r : this->relations) {
        if ((r->name == rel_id) ||
            (!r->alias.empty() && (r->alias == rel_id)))
            return r;
    }
    if (where && (whereClause != nullptr)) { // Check if relation exists in WHERE clause
        if (auto r = whereClause->getRelation(rel_id)) {
            return std::dynamic_pointer_cast<RelationExpr>(r);
        }
    }
    return nullptr;
}

// Returns true if `this` is a sub-query that appears in a WHERE clause
bool QueryGraph::appearsInWhereClause() const {
    if (this->parent == nullptr)
        return false;
    auto where = parent->whereClause;
    if (where == nullptr)
        return false;
    // Get query name as used in parent
    std::string name = "#Q_" + std::to_string(this->id);
    return where->contains(name);
}

// Checks if `this` query returns the given attribute
bool QueryGraph::outputs(std::shared_ptr<AttributeExpr> a) {
    for (auto r : this->projections) {
        if (r->alias.empty()) {
            if (r->hasType(Expression::AttributeRef)) {
                auto att = std::dynamic_pointer_cast<AttributeExpr>(r);
                bool flag = (att->name == a->name);
                if (!a->relation.empty()) {
                    // TODO (john): Fix
//                    assert(!att->relation.empty()); // Make sure schema is resolved
//                    return flag && (att->relation == a->relation);
                    return true;
                }
                else if (flag)
                    return true;
            }
        }
        else if (r->alias == a->name) {
            return true;
        }
    }
    return false;
}

bool QueryGraph::containsStar() const {
    for (auto a : this->projections) {
        if (a->hasType(Expression::Star)) {
            return true;
        }
    }
    return false;
}

bool QueryGraph::contains(const std::string &rel_name) const {
    for (auto r : this->relations) {
        if (r->name == rel_name || r->alias == rel_name) {
            return true;
        }
    }
    return false;
}

void QueryGraph::printHierarchy(int indent) const {
    if (get_rank()==0) {
        std::string prefix = indent==0 ? "" : "-";
        for (int i=0; i<indent; i++)
            std::cout << "  ";
        std::cout << prefix << "QUERY " << this->id << std::endl;
        for (auto c : this->children) {
            c->printHierarchy(indent+1);
        }
    }
}

void QueryGraph::printInfo() const {
    if (get_rank()==0) {
        std::cout << "===========================" << std::endl;
        std::cout << "QUERY GRAPH INFO (ID: " << this->id << ")" << std::endl;
        std::cout << "===========================" << std::endl;
        std::string proj = this->distinct ? "PROJECTIONS (DISTINCT):" : "PROJECTIONS:";
        std::cout << proj << std::endl;
        for (auto p : this->projections)
            std::cout << p->toString() << std::endl;
        std::cout << "===========================" << std::endl;
        std::cout << "RELATIONS: " << std::endl;
        for (auto r : this->relations) {
            std::cout << r->name;
            if (!r->alias.empty())
                std::cout << " (alias: " << r->alias << ")";
            std::cout << std::endl;
        }
        std::cout << "===========================" << std::endl;
        std::cout << "JOINS: " << std::endl;
        for (auto j : this->joins) {
            j->print();
        }
        std::cout << "===========================" << std::endl;
        std::cout << "FILTERS: " << std::endl;
        for (auto f : this->filters)
            f->print();
        std::cout << "===========================" << std::endl;
        if (this->groupBy != nullptr) {
            std::cout << "GROUP-BY: " << std::endl;
            std::cout << "Keys: ";
            for (auto key: this->groupBy->keys)
                std::cout << key->toString() << "   ";
            std::cout << std::endl;
            std::cout << "Aggregation expressions: ";
            for (auto att: this->groupBy->aggs)
                std::cout << att->toString() << "   ";
            std::cout << std::endl;
            std::cout << "Having expressions: ";
            if (this->groupBy->having != nullptr)
                std::cout << this->groupBy->having->toString() << "   ";
            std::cout << std::endl;
            std::cout << "===========================" << std::endl;
        }
        if (this->orderBy != nullptr) {
            std::cout << "ORDER-BY: " << std::endl;
            std::cout << "Keys: ";
            for (auto key: this->orderBy->keys)
                std::cout << key->toString() << "   ";
            std::cout << std::endl;
            std::cout << "Direction: ";
            std::string dir;
            for (auto key: this->orderBy->asc)
                dir = (key == true) ? "ASC" : "DESC";
            std::cout << dir << "   ";
            std::cout << std::endl;
            std::cout << "===========================" << std::endl;
        }
        for (auto c : this->children) {
            c->printInfo();
        }
    }
}

void QueryGraph::collectQueries(std::map<int, std::shared_ptr<QueryGraph>> &m) {
    m.insert(std::make_pair(this->id, shared_from_this()));
    for (auto q : this->children)
        q->collectQueries(m);
    if (this->sibling != nullptr)
        this->sibling->collectQueries(m);
}

std::shared_ptr<FilterExpr> QueryGraph::isFilter(const std::shared_ptr<Expression> expr) const {
    if (expr->isComparison()) {
        std::vector <std::shared_ptr<Expression>> attributes;
        expr->collectAttributes(attributes);
        std::set <std::string> s;
        for (auto a: attributes) {
            assert(!a->relation.empty());
            s.insert(a->relation);
        }
        if (s.size() == 1) {
            return std::make_shared<FilterExpr>(FilterExpr{expr->expressionType, attributes, expr});
        }
    }
    return nullptr;
}

std::shared_ptr<JoinExpr>  QueryGraph::isJoin(std::shared_ptr<Expression> expr, JoinExpr::JoinType joinType) const {
    if (expr->hasType(Expression::In)) { // Semi-join
        assert(expr->expressions.size() == 2);
        assert(expr->expressions[0]->hasType(Expression::AttributeRef) &&
                expr->expressions[1]->hasType(Expression::RelationRef));
        auto left = this->getRelation(expr->expressions[0]->relation, false);
        auto right = std::dynamic_pointer_cast<RelationExpr>(expr->expressions[1]);
        return std::make_shared<JoinExpr>(JoinExpr{JoinExpr::Semi, expr, left, right});
    }
    else {  // Extract relations
        assert(expr->isComparison());
        std::vector<std::shared_ptr<Expression>> attributes;
        std::set<std::string> relations;
        expr->collectAttributes(attributes);
        for (auto a : attributes) {
            assert(!a->relation.empty());
            relations.insert(a->relation);
        }
        if (relations.size()==2) { // We have a join
            std::vector<std::string> t(relations.begin(), relations.end());
            std::shared_ptr<RelationExpr> left, right;
            for (auto a : attributes) {
                if (a->relation == t[0]) {
                    left = this->getRelation(t[0]);
                }
                else {
                    assert(a->relation == t[1]);
                    right = this->getRelation(t[1]);
                }
            }
            assert(left != nullptr);
            assert(right != nullptr);
            return std::make_shared<JoinExpr>(JoinExpr{joinType, expr, left, right});
        }
        else if (relations.size()>2) {
            leaderLogging("[ERROR] Unsupported join predicate: " + expr->toString() + "\n.", true);
        }
    }
    return nullptr;
}

int QueryGraph::getRelId(const std::string &name_or_alias) const {
    for (int i=0; i<this->relations.size(); i++)
        if (this->relations[i]->identifiedBy(name_or_alias))
            return i;
    leaderLogging("[ERROR] Relation " + name_or_alias + " not found\n", true);
    return -1;
}

void QueryGraph::collectFiltersAndJoins(std::shared_ptr<Expression> expr, JoinExpr::JoinType joinType) {
    std::vector<std::shared_ptr<Expression>> predicates;
    expr->flatten(predicates);
    for (auto p : predicates) {
        if (auto f = this->isFilter(p)) {
            this->filters.push_back(f);
            int rId = this->getRelId(f->attributes[0]->relation);
            f->rId = rId;
        }
        else if (auto j = this->isJoin(p, joinType)) {
            this->joins.push_back(j);
            std::string s1 = j->left->alias.empty() ? j->left->name : j->left->alias;
            std::string s2 = j->right->alias.empty() ? j->right->name : j->right->alias;
            j->rId1 = this->getRelId(s1);
            j->rId2 = this->getRelId(s2);
        }
    }
}

void QueryGraph::collectTableFunctions() {
    if (groupBy != nullptr)
        return;
    for (auto e : projections) {
        if (e->hasType(Expression::FunctionRef)) {
            auto tf = std::make_shared<TableFunctionExpr>(TableFunctionExpr(e));
            tFunctions.push_back(tf);
        }
    }
}

unsigned QueryGraph::numOp() const {
    return filters.size() +
            joins.size() +
            (unsigned) distinct +
            (groupBy != nullptr) +
            (orderBy != nullptr) +
            tFunctions.size();
}

unsigned QueryGraph::getDistinctId() const {
    return filters.size() +
           joins.size() +
           (groupBy != nullptr) +
           tFunctions.size();
}

void QueryGraph::collectFiltersAndJoins() {
    // Check predicates in WHERE clause
    if (this->whereClause != nullptr) {
        this->collectFiltersAndJoins(this->whereClause);
    }
    // Check predicates in FROM clause
    if (this->fromClause != nullptr) {
        this->collectFiltersAndJoins(this->fromClause, this->jType);
        // Merge predicates in FROM and WHERE
        if (whereClause) {
            auto root = std::make_shared<Expression>(Expression::And);
            root->expressions.push_back(whereClause);
            root->expressions.push_back(fromClause);
            this->whereClause->parent = root;
            this->fromClause->parent = root;
            this->whereClause = root;
        }
        else {
            whereClause = fromClause;
        }
    }
    // Check children and sibling
    for (auto c : this->children)
        c->collectFiltersAndJoins();
    if (this->sibling != nullptr)
        this->sibling->collectFiltersAndJoins();
}

bool QueryGraph::hasAggregation() const {
    for (auto p : this->projections) {
        if (p->hasType(Expression::FunctionRef))
            return true;
    }
    return false;
}

// Construct the edges
void QueryGraph::constructEdges(const Database& db)
{
    if (this->whereClause) {
        std::vector<std::shared_ptr<Expression>> tmpRels;
        this->whereClause->collectTempRelations(tmpRels);
        for (auto r : tmpRels) {
            auto rel = std::dynamic_pointer_cast<RelationExpr>(r);
            this->relations.push_back(rel);
        }
    }
    collectFiltersAndJoins();
    collectTableFunctions();
    setOpIds();
    setBaseCardinalities(db);
}

void QueryGraph::setOpIds() {
    for (auto c : children)
        c->setOpIds();
    if (sibling != nullptr)
        sibling->setOpIds();
    // Step 0: Filter ids
    for (int i=0; i<filters.size(); i++) {
        filters[i]->id = i;
        filters[i]->selPredicate->opId = i;
    }
    int n = this->filters.size();
    // Step 1: Join ids
    for (int i=0; i<joins.size(); i++) {
        joins[i]->id = i+n;
        joins[i]->predicate->opId = i+n;
    }
    n += joins.size();
    // Step 2: Group-by id
    if (groupBy) {
        groupBy->id = n;
        for (auto a : groupBy->aggs) {
            if (!a->hasType(Expression::AttributeRef)) {
                a->opId = n;
            }
        }
        if (groupBy->having) {
            groupBy->having->opId = n;
        }
        n++;
    }
    // Step 3: Order-by id
    if (orderBy) {
        orderBy->id = n++;
    }
    // Step 4: Table function ids
    for (auto f : tFunctions) {
        f->id = n;
        f->exp->opId = n++;
    }
}

void QueryGraph::setBaseCardinalities(const Database &db) {
    for (auto q : children)
        q->setBaseCardinalities(db);
    if (sibling != nullptr)
        sibling->setBaseCardinalities(db);
    for (auto r : relations) {
        if (!r->subQuery) {
            auto rel = db.getRelation(r->name);
            r->cardinality = rel->cardinality();
        }
    }
};

// Prints information about the graph
void QueryGraph::print() const {
    if (get_rank()==0) {
        if (!this->children.empty()) { // Print hierarchy of sub-queries
            std::cout << "===========================" << std::endl;
            std::cout << "QUERY HIERARCHY" << std::endl;
            std::cout << "===========================" << std::endl;
            this->printHierarchy();
            std::cout << "===========================" << std::endl;
        }
        // Print query info
        this->printInfo();
        // Print next query
        if (this->sibling != nullptr)
            this->sibling->printInfo();
    }
}