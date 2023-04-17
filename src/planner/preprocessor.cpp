#include "../../include/planner/preprocessor.h"

/***
 * Constructor
 */
Preprocessor::Preprocessor() { }

/***
 * Destructor
 */
Preprocessor::~Preprocessor() { }

// Checks if state is empty
bool Preprocessor::stateIsEmpty() const {
    return this->qc==0 &&
            this->queue.empty() &&
            this->qIndex.empty() &&
            this->sIndex.empty() &&
            this->a2id.empty();
}

// Sets the name of the temporary relation that represents the output of the given sub-query
void Preprocessor::setAlias(Parser::SQLQuery &subQuery) {
    auto p = this->sIndex.find(subQuery.tree);
    if (p != this->sIndex.end()) {
        auto v = p->second;
        for (auto r : v) {
            if (r->alias.empty())
                r->name += std::to_string(subQuery.id);
            else {  // Update alias-to-id mapping
                auto p = this->a2id.insert(std::make_pair(r->alias, subQuery.id));
                assert(p.second);
            }
        }
    }
//    else {
//        leaderLogging("[ERROR] Sub-query with id " + std::to_string(subQuery.id) + " not found.\n", true);
//    }
}

// Retrieves the graph of the query with given `name`
std::shared_ptr<QueryGraph> Preprocessor::resolveQuery(std::shared_ptr<RelationExpr> rel) const {
    // The given name may be and alias or a temporary relation name of the form 'TMP_PREFIXid'
    std::string name = rel->alias.empty() ? rel->name : rel->alias;
    auto p = this->a2id.find(name);
    int qId = (p != this->a2id.end()) ? qId = p->second :
                                        std::stoi(name.substr(name.find(TMP_PREFIX) + TMP_PREFIX.size()));
    auto q = this->qIndex.find(qId);
    if (q != this->qIndex.end()) {
        return q->second;
    }
    std::string msg = "[ERROR]: Query with id " + std::to_string(qId) + " not found.\n";
    leaderLogging(msg, true);
    return nullptr;
}

// Resolves the schema of the query result
void Preprocessor::resolveSchema(std::shared_ptr<QueryGraph> query, Database &db,
                                 std::shared_ptr<RelationExpr> relation) const {
    if (!relation && !query->containsStar())
        return;
    // Resolve star
    std::vector<std::shared_ptr<Expression>> attributes;
    for (auto r : query->relations) {
        if (r->subQuery) { // Sub-query
            std::shared_ptr<QueryGraph> child = this->resolveQuery(r);
            this->resolveSchema(child, db, r);
            for (auto a : child->projections) {
                a->qID = child->id; // Keep track of the query origin
                attributes.push_back(a);
            }
        }
        else { // Base relation
            auto rel = db.getRelation(r->name);
            assert(rel != nullptr);
            auto rel_name = r->alias.empty() ? r->name : r->alias;
            auto cols = rel->getCols(true);
            for (auto c : cols) { // Create new attribute expression
                auto att = std::make_shared<AttributeExpr>(AttributeExpr(c.c_str(),
                                                                                                rel_name.c_str()));
                attributes.push_back(att);
            }
        }
    }
    // Update projections
    std::vector<std::shared_ptr<Expression>> new_proj;
    for (auto p : query->projections) {
        if (p->hasType(Expression::Star)) { // Replace star with list of attributes
            for (auto a : attributes)
                new_proj.push_back(a);
        }
        else { // Keep expression as is
            new_proj.push_back(p);
        }
    }
    query->projections = new_proj;
    if ((relation != nullptr) && relation->subQuery) {
        if (relation->customSchema) {  // Sub-query with user-defined schema
            assert(query->projections.size() == relation->schema.size());
            for (int i = 0; i < query->projections.size(); i++) {
                assert(relation->schema[i]->alias.empty());
                query->projections[i]->alias = relation->schema[i]->name;
            }
        }
        else if (relation->schema.empty()) {
            relation->schema = query->projections;
        }
        else {
            assert(relation->schema.size() == query->projections.size());
        }
    }
}

// Checks if the given attribute exists in a base relation or the result of a sub-query
bool Preprocessor::validate(std::shared_ptr<AttributeExpr> att, std::shared_ptr<RelationExpr> rel,
                            Database &db) const {
    if (rel->subQuery) { // Sub-query
        std::shared_ptr<QueryGraph> query = this->resolveQuery(rel);
        this->resolveSchema(query, db, rel);  // Make sure the schema is resolved
        if (query->outputs(att))
            return true;
    }
    else if (db.contains(att->name, rel->name)) { // Base relation
            return true;
    }
    return false;
}

// Checks if a relation with name `rName` exists in the database
void Preprocessor::validate(std::string rName, Database &db) const {
    if (!db.containsRelation(rName)) {
        std::string msg = "[ERROR] Relation '" + rName + "' does not exist in the database.\n";
        leaderLogging(msg, true);
    }
}

// Checks if all base relations in the query exist in the database
void Preprocessor::validateRelations(std::shared_ptr<QueryGraph> query, Database &db) const {
    leaderLogging("[INFO] Validating relations...\n");
    std::set<std::string> n2a;
    std::map<std::string, std::string> a2n;
    for (auto r : query->relations) {
        auto p1 = n2a.insert(r->name + r->alias);
        if (!p1.second) { // 'name+alias' exists
            std::string msg = "[ERROR] Relation '" + r->name + "' appears more than once in the query.\n";
            leaderLogging(msg, true);
        }
        if (!r->alias.empty()) {
            auto p2 = a2n.insert(std::make_pair(r->alias, r->name));
            if (!p2.second) { // Alias exists
                std::string msg = "[ERROR] Alias '" + r->alias + "' appears more than once in the query.\n";
                leaderLogging(msg, true);
            }
        }
        if (!r->subQuery) { // Check if base relation exists in the database
            this->validate(r->name, db);
            // Set schema
            auto br = db.getRelation(r->name);
            auto cols = br->getCols(true);
            for (auto a : cols) {
                auto n = r->alias.empty() ? r->name : r->alias;
                auto att = std::make_shared<AttributeExpr>(AttributeExpr(a.c_str(),
                                                                                                n.c_str()));
                r->schema.push_back(att);
            }
        }
    }
}

// Check if the given attribute expression is valid w.r.t. the given query and database
void Preprocessor::validate(std::shared_ptr<AttributeExpr> a,
                            std::shared_ptr<QueryGraph> query,
                            Database &db, bool sortingKey) const {
    if (a->relation.empty()) { // Orphan attribute
        // Get candidate relations from database
        std::vector<std::shared_ptr<Relation>> from_db;
        db.getRelations(a->name, from_db);
        std::vector<std::shared_ptr<RelationExpr>> candidates;
        // Filter out those relations that do not appear in the query
        for (auto r : query->relations) {
            for (auto rb : from_db) {
                if (r->name == rb->getName())
                    candidates.push_back(r);
            }
        }
        auto m = candidates.size();   // Number of candidate base relations
        if (m == 1) { // Attribute exists in a base relation
            a->relation = candidates[0]->alias.empty() ? candidates[0]->name : candidates[0]->alias;
            return;
        }
        else if (m > 1){ // Attribute is ambiguous
            leaderLogging("[ERROR] Attribute '" + a->name + "' is ambiguous.\n", true);
        }
        else { // Attribute may appear in the result of a sub-query
            int cnt = 0;
            for (auto rel : query->relations) { // For all relation in the FROM clause
                if (this->validate(a, rel, db)) {
                    a->relation = rel->alias.empty() ? rel->name : rel->alias;
                    cnt++;
                }
            }
            if (cnt==1)
                return;
            else if (cnt>1)
                leaderLogging("[ERROR] Attribute '" + a->name + "' is ambiguous.\n", true);
        }
    }
    else {// Attribute has an associated relation name or alias
        std::vector<std::shared_ptr<RelationExpr> > v;
        query->getRelations(a->relation, v);
        for (auto rel : v) {
            if (this->validate(a, rel, db)) // Relation exists in FROM clause
                return;
        }
    }
    // Attribute is not resolved, check if it is defined in the SELECT clause or belongs to a relation in a parent query
    if (sortingKey) { // This is a sorting key
        if (query->outputs(a)) {
            a->relation = TMP_PREFIX + "0";  // Attribute is defined in the SELECT clause
            return;
        }
    }
    else if (query->appearsInWhereClause()) { // This is a sub-query that appears in a WHERE clause
        this->validate(a, query->parent, db);  // Check if attribute's relation appears in parent
        return;
    }
    leaderLogging("[ERROR] Attribute '" + a->name + "' could not be resolved.\n", true);
}

// Validates all attributes in the given expression
void Preprocessor::validate(std::shared_ptr<Expression> &expr,
                            std::shared_ptr<QueryGraph> &query,
                            Database &db, bool sortingKey) const {
    assert(expr != nullptr);
    if (!expr->hasType(Expression::AttributeRef)) {
        for (auto e: expr->expressions)
            this->validate(e, query, db);
    }
    else { // Expression is an attribute
        auto a = std::dynamic_pointer_cast<AttributeExpr>(expr);
        this->validate(a, query, db, sortingKey);
    }
}

// Validates attribute names and updates their relation information accordingly
void Preprocessor::validateAttributes(std::shared_ptr<QueryGraph> query, Database &db) const {
    leaderLogging("[INFO] Validating attributes...\n");
    // Check attributes in SELECT clause
    for (auto j: query->projections) {
        this->validate(j, query, db);
    }
    // Check attributes in FROM clause
    if (query->fromClause != nullptr) {
        this->validate(query->fromClause, query, db);
    }
    // Check attributes in WHERE clause
    if (query->whereClause != nullptr) {
        this->validate(query->whereClause, query, db);
    }
    // Check attributes in GROUP-BY clause
    if (query->groupBy != nullptr) {
        for (auto k: query->groupBy->keys) {
            this->validate(k, query, db);
        }
        for (auto a: query->groupBy->aggs) {
            this->validate(a, query, db);
        }
        auto having = query->groupBy->having;
        if (having != nullptr) {
            this->validate(having, query, db);
        }
    }
    // Check attributes in ORDER-BY clause
    if (query->orderBy != nullptr) {
        for (auto k: query->orderBy->keys) {
            this->validate(k, query, db, true);
        }
    }
}

// Validates IN operations in the given expression
void Preprocessor::validateSemiJoins(std::shared_ptr<Expression> expr, std::shared_ptr<QueryGraph> query,
                                     Database &db) const {
    for (auto e : expr->expressions) {
        this->validateSemiJoins(e, query, db);
    }
    if (expr->hasType(Expression::In)) { // Make sure the two operands have the same schema
        auto left = expr->expressions[0];
        auto right = expr->expressions[1];
        assert(left->hasType(Expression::AttributeRef));
        assert(right->hasType(Expression::RelationRef));
        auto r = std::dynamic_pointer_cast<RelationExpr>(right);
        assert(r->subQuery);
        std::shared_ptr<QueryGraph> child = this->resolveQuery(r);
        this->resolveSchema(child, db, r);
        if (child->projections.size() == 2) // Omit SEL column
            return;
        leaderLogging("[ERROR] Incompatible semi-join operands.\n", true);
        // TODO: Check plaintext data types
    }
}

// Validates IN operations in the query
void Preprocessor::validateSemiJoins(std::shared_ptr<QueryGraph> query, Database &db) const {
    leaderLogging("[INFO] Validating semi-joins...\n");
    if (query->whereClause != nullptr)
        this->validateSemiJoins(query->whereClause, query, db);
}

// Validates query and updates incomplete attribute expressions accordingly
void Preprocessor::validate(std::shared_ptr<QueryGraph> query, Database &db) {
    std::string msg = "[INFO] Validating query with ID " + std::to_string(query->id) + "...\n";
    leaderLogging(msg);
    this->validateRelations(query, db);
    this->validateAttributes(query, db);
    this->validateSemiJoins(query, db);
    query->validated = true;
}

// Processes the given SQL expression
std::shared_ptr<Expression> Preprocessor::processExpression(const hsql::Expr *expr, Parser::SQLQuery& query,
                                                            std::shared_ptr<Expression> parent) {
    assert(expr != nullptr);
    Expression::ExpressionType eType;
    // Identify the root of the given AST
    switch (expr->opType) {
        case hsql::OperatorType::kOpPlus:
            eType = Expression::Addition;
            break;
        case hsql::OperatorType::kOpMinus:
            eType = Expression::Subtraction;
            break;
        case hsql::OperatorType::kOpAsterisk:
            eType = Expression::Multiplication;
            break;
        case hsql::OperatorType::kOpSlash:
            eType = Expression::Division;
            break;
        case hsql::OperatorType::kOpEquals:
            eType = Expression::Equal;
            break;
        case hsql::OperatorType::kOpNotEquals:
            eType = Expression::NotEqual;
            break;
        case hsql::OperatorType::kOpLess:
            eType = Expression::LessThan;
            break;
        case hsql::OperatorType::kOpLessEq:
            eType = Expression::LessThanEqual;
            break;
        case hsql::OperatorType::kOpGreater:
            eType = Expression::GreaterThan;
            break;
        case hsql::OperatorType::kOpGreaterEq:
            eType = Expression::GreaterThanEqual;
            break;
        case hsql::OperatorType::kOpAnd:
            eType = Expression::And;
            break;
        case hsql::OperatorType::kOpOr:
            eType = Expression::Or;
            break;
        case hsql::OperatorType::kOpIn:
            eType = Expression::In;
            break;
        case hsql::OperatorType::kOpExists:
            eType = Expression::Exists;
            break;
        default:
            eType = Expression::Unknown;
    }

    if (eType == Expression::Unknown) { // Expression is not an operation
        switch (expr->type) {
            case hsql::kExprColumnRef: // Expression is a column name
            {
                std::string table = expr->table == nullptr ? "" : expr->table;
                return std::make_shared<AttributeExpr>(AttributeExpr(expr->name,
                                                                        expr->table,
                                                                        expr->alias,
                                                                        parent));
            }
            case hsql::kExprLiteralString: // Expression is a string
                return std::make_shared<ConstantExpr>(ConstantExpr(expr->name, expr->alias, parent));
            case hsql::kExprLiteralInt: // Expression is an integer
                return std::make_shared<ConstantExpr>(ConstantExpr(expr->ival, expr->alias, parent));
            case hsql::kExprLiteralFloat: // Expression is a float
                return std::make_shared<ConstantExpr>(ConstantExpr(expr->fval, expr->alias, parent));
            case hsql::ExprType::kExprStar: // Expression is a '*' (all)
                return std::make_shared<Expression>(Expression(Expression::Star, expr->alias, parent));
            case hsql::kExprFunctionRef: // Expression is a function
            {
                auto func = std::make_shared<FunctionExpr>(FunctionExpr(expr->name,
                                                                        expr->alias,
                                                                        expr->distinct,
                                                                        parent));
                // Collect input of the aggregation function
                func->expressions.push_back(processExpression(expr->exprList->at(0), query, func));
                return func;
            }
            case hsql::kExprSelect: // Expression is a sub-query
            {
                // Anonymous sub-queries are assigned an id of the form "#Q_"
                query.subTrees.push_back(expr->select); // Update list of sub-queries
                auto tRef = std::make_shared<RelationExpr>(RelationExpr(TMP_PREFIX.c_str(),
                                                                                                expr->alias,
                                                                                                true,
                                                                                                nullptr,
                                                                                                parent));

                // Cache the sub-query to update its id afterwards
                std::vector<std::shared_ptr<RelationExpr>> v = {tRef};
                this->sIndex.insert(std::make_pair(expr->select, v));
                return tRef;
            }
            default:
                leaderLogging("[ERROR] Unsupported SQL expression.\n", true);
        }
    }
    // Initialize expression
    auto expression = std::make_shared<Expression>(Expression(eType));
    expression->parent = parent;
    // Collect sub-expressions (operands)
    if (expr->expr != nullptr) { // First operand
        expression->addExpression(processExpression(expr->expr, query, expression));
    }
    if (expr->expr2 != nullptr) { // First operand
        expression->addExpression(processExpression(expr->expr2, query, expression));
    }
    if (expr->exprList != nullptr) { // N-ary operation or sub-query with alias
        for (auto exp : expr->exprList[0]) {
            auto child = processExpression(exp, query, expression);
            expression->addExpression(child);
        }
    }
    if (expr->select != nullptr) { // Sub-query without alias
        assert(expr->alias == nullptr);
        // Anonymous sub-queries are assigned an id of the form "#Q_"
        query.subTrees.push_back(expr->select); // Update list of sub-queries
        auto child = std::make_shared<RelationExpr>(RelationExpr(TMP_PREFIX.c_str(),
                                                                                        nullptr,
                                                                                        true,
                                                                                        nullptr,
                                                                                        expression));
        expression->addExpression(child);
        // Cache the sub-query to update its id afterwards
        std::vector<std::shared_ptr<RelationExpr>> v = {child};
        this->sIndex.insert(std::make_pair(expr->select, v));
    }
    // Make sure the expression has been parsed correctly
    assert(eType != Expression::Unknown);
    return expression;
}

// Processes a reference base relation
std::shared_ptr<RelationExpr> Preprocessor::processTableRef(const hsql::TableRef* tableRef) {
    char* alias = tableRef->alias != nullptr ? tableRef->alias->name : nullptr;
    return std::make_shared<RelationExpr>(RelationExpr(tableRef->name, alias));
}

// Processes a join reference
std::shared_ptr<JoinExpr> Preprocessor::processJoin(const hsql::JoinDefinition* joinRef, Parser::SQLQuery &query) {
    // Extract join type
    JoinExpr::JoinType jType;
    switch (joinRef->type) {
        case hsql::kJoinInner:
            jType = JoinExpr::Inner;
            break;
        case hsql::kJoinFull:
            jType = JoinExpr::Full;
            break;
        case hsql::kJoinLeft:
            jType = JoinExpr::LeftOuter;
            break;
        case hsql::kJoinRight:
            jType = JoinExpr::RightOuter;
            break;
        case hsql::kJoinNatural:
        case hsql::kJoinCross: {
            leaderLogging("[ERROR] Unsupported join type.\n", true);
        }
    }
    // Process join predicate
    auto pred = processExpression(joinRef->condition, query);
    // Extract left and right input relations
    assert((joinRef->left->name != nullptr) && (joinRef->right->name != nullptr));
    auto left = processTableRef(joinRef->left);
    auto right = processTableRef(joinRef->right);
    return std::make_shared<JoinExpr>(JoinExpr {jType, pred, left, right});
}

// Processes SELECT clause
void Preprocessor::processSelectClause(Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc) {
    auto selectStatement = query.tree;
    if (selectStatement->selectList != nullptr) {
        queryDesc->distinct = selectStatement->selectDistinct;
        for (auto expr: *selectStatement->selectList)
            queryDesc->projections.push_back(processExpression(expr, query));
    }
}

// Processes FROM clause
void Preprocessor::processFromClause(Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc) {
    auto tableRef = query.tree->fromTable;
    // Collect all base tables, sub-queries, and join descriptions
    if (tableRef->list != nullptr) {  // Two or more tables
        for (auto t : *tableRef->list) {
            // TODO (john): fix
            assert(t->name != nullptr);
            queryDesc->relations.push_back(processTableRef(t));
        }
    }
    else if (tableRef->name != nullptr){ // Only one table
        queryDesc->relations.push_back(processTableRef(tableRef));
    }
    else if (tableRef->select != nullptr) { // Sub-query
        query.subTrees.push_back(tableRef->select);
        // Anonymous sub-queries are assigned an id of the form "#Q_"
        std::string alias = tableRef->alias != nullptr ? tableRef->alias->name : TMP_PREFIX;
        std::vector<char*>* schema = tableRef->alias != nullptr ? tableRef->alias->columns : nullptr;
        auto tRef = std::make_shared<RelationExpr>(RelationExpr(alias.c_str(),
                                                                                alias.c_str(),
                                                                            true, schema));
        queryDesc->relations.push_back(tRef);
        // Cache the sub-query to update its id afterwards
        std::vector<std::shared_ptr<RelationExpr>> v = {tRef};
        this->sIndex.insert(std::make_pair(tableRef->select, v));
    }
    else if (tableRef->join != nullptr) { // FROM clause with join description
        auto jDesc = processJoin(tableRef->join, query);
        queryDesc->fromClause = jDesc->predicate;
        queryDesc->jType = jDesc->jType;
//        queryDesc->joins.push_back(jDesc);
        // Update the list of relations
        queryDesc->relations.push_back(jDesc->left);
        queryDesc->relations.push_back(jDesc->right);
    }
    else {
        leaderLogging("[ERROR] Unsupported FROM clause.\n");
        exit(-1);
    }
}

// Processes WHERE clause
void Preprocessor::processWhereClause(Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc) {
    if (query.whereClause() != nullptr) { // If there is a WHERE clause
        queryDesc->whereClause = processExpression(query.whereClause(), query);
    }
}

// Processes GROUP BY clause
void Preprocessor::processGroupByClause(Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc) {
    auto selectStatement = query.tree;
    if (selectStatement->groupBy != nullptr) {
        queryDesc->groupBy = std::make_shared<GroupByExpr>(GroupByExpr {});
        // Collect grouping keys
        for (int i = 0; i < selectStatement->groupBy[0].columns[0].size(); ++i) {
            auto column = selectStatement->groupBy[0].columns[0][i];
            queryDesc->groupBy->keys.push_back(std::make_shared<AttributeExpr>(AttributeExpr(column->name,
                                                                                            column->table,
                                                                                            column->alias)));
        }
        // Collect HAVING information
        if (selectStatement->groupBy->having != nullptr) {
            queryDesc->groupBy->having = processExpression(selectStatement->groupBy->having,query);
        }
        // Collect aggregation expressions in SELECT clause
        if (selectStatement->selectList != nullptr) {
            for (auto expr : *selectStatement->selectList) {
                if (expr->type == hsql::kExprFunctionRef) {
                    std::shared_ptr<Expression> e = processExpression(expr, query);
                    assert(e->expressions.size()>0);
                    queryDesc->groupBy->aggs.push_back(e);
                }
            }
        }
    }
}

// Processes ORDER BY clause
void Preprocessor::processOrderByClause(const Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc) {
    auto selectStatement = query.tree;
    if (selectStatement->order != nullptr) {
        queryDesc->orderBy = std::make_shared<OrderByExpr>(OrderByExpr {});
        for (int i = 0; i < selectStatement->order[0].size(); ++i) {
            hsql::OrderDescription *oDesc = selectStatement->order[0][i];
            queryDesc->orderBy->asc.push_back(oDesc->type == hsql::OrderType::kOrderAsc ? true : false);
            queryDesc->orderBy->keys.push_back(std::make_shared<AttributeExpr>(AttributeExpr(oDesc->expr->name,
                                                                                            oDesc->expr->table,
                                                                                            oDesc->expr->alias)));
        }
    }
    // Collect LIMIT information
    // TODO: LIMIT may exist without an ORDER-BY clause
    if (selectStatement->limit != nullptr) {
        queryDesc->orderBy->limit = selectStatement->limit->limit->ival;
    }
}
// Processes ORDER BY clause
void Preprocessor::processWithClause(const Parser::SQLQuery &query, std::shared_ptr<QueryGraph> queryDesc) {
    if (query.withClause() != nullptr) {
        std::cout << "WITH SIZE: " << query.withClause()->size() << std::endl;
    }
}


// Collects information about the given query
std::shared_ptr<QueryGraph> Preprocessor::processQuery(Parser::SQLQuery& query, int id,
                                                       std::shared_ptr<QueryGraph> parent) {
    query.id = id;
    auto qd = std::make_shared<QueryGraph>(QueryGraph(id, parent));
    leaderLogging("[INFO] Processing SELECT clause.\n");
    processSelectClause(query, qd);
    leaderLogging("[INFO] Processing FROM clause.\n");
    processFromClause(query, qd);
    leaderLogging("[INFO] Processing WHERE clause.\n");
    processWhereClause(query, qd);
    leaderLogging("[INFO] Processing GROUP BY clause.\n");
    processGroupByClause(query, qd);
    leaderLogging("[INFO] Processing ORDER BY clause.\n");
    processOrderByClause(query, qd);
    // Add sub-queries to the queue
    for (auto s : query.subTrees) {
        Parser::SQLQuery sq = {.tree = s};
        this->queue.push(std::make_pair(sq, qd));
    }
    // Update query index
    auto p = this->qIndex.insert(std::make_pair(id, qd));
    assert(p.second != false);
    // Set query alias
    this->setAlias(query);
    auto a=sIndex.find(query.tree);
    if (a!=sIndex.end()) {
        qd->relAlias = a->second;
    }
    return qd;
}

// Processes all queries in the queue
void Preprocessor::processQueue() {
    while (!this->queue.empty()) {
        // Collect information about sub-query
        auto next = this->queue.front();
        Parser::SQLQuery q = next.first;
        auto p = next.second;
        auto qd = this->processQuery(q, this->qc++, p);
        // Update parent
        if (p != nullptr)
            p->children.push_back(qd);
        // Update queue
        this->queue.pop();
    }
}

/***
 * Collects information from the parsed SQL query.
 * @param query - The parsed SQL query.
 * @return The query graph.
 */
std::shared_ptr<QueryGraph> Preprocessor::collectInfo(Parser::SQLQuery &query) {
    assert(stateIsEmpty());
    // Collect information about the query
    auto root = this->processQuery(query, this->qc++);
    // Process sub-queries
    this->processQueue();
    // Process queries in WITH clause
    if (query.withClause() != nullptr) {
        for (auto q : *query.withClause()) {
            assert(q->alias != nullptr);
            // Find corresponding relation and update sub-query index
            std::vector<std::shared_ptr<RelationExpr>> v;
            root->getRelations(q->alias, v, true);
            for (auto r : v)
                r->subQuery = true;
            this->sIndex.insert(std::make_pair(q->select, v));
            Parser::SQLQuery sq = {.tree = q->select};
            auto next = this->processQuery(sq, this->qc++);
            next->relAlias = v;
            root->sibling = next;
            this->processQueue();
        }
    }
    return root;
}

// Applies semantic analysis and transforms the given query into an equivalent one
void Preprocessor::preprocess(std::shared_ptr<QueryGraph> query, Database &db) {
    // TODO: More transformations are described by R. A. Ganski et al. in "Optimization for Nested SQL Queries
    //  Revisited", SIGMOD Rec., 1987.

    // If there is an EXISTS clause with a sub-query, let q1, whose WHERE clause contains an attribute of a relation
    // in the outer query, let q0, merge q1's FROM and WHERE clauses into q0's FROM and WHERE clauses and discard q1
    if (query->appearsInWhereClause()) {
        // Get temporary relation name as used in parent
        std::string name = TMP_PREFIX + std::to_string(query->id);
        // Get parent's expression that corresponds to the sub-query
        auto r = query->parent->whereClause->getRelation(name);
        assert(r != nullptr);
        if (r->parent->hasType(Expression::Exists)) {
            if (!query->hasAggregation()) {
                // Get attributes in WHERE clause of the sub-query
                auto w = query->whereClause;
                std::vector<std::shared_ptr<Expression>> va;
                w->collectAttributes(va);
                for (auto att: va) {
                    auto a = std::dynamic_pointer_cast<AttributeExpr>(att);
                    // If there is an attribute that belongs to the parent relation
                    if (!query->contains(a->relation) && query->parent->contains(a->relation)) {
                        // Merge child's and parent's FROM clauses
                        for (auto r: query->relations)
                            query->parent->relations.push_back(r);
                        // Merge child's and parent's WHERE clauses
                        auto grandparent = r->parent->parent;
                        grandparent->expressions.erase(std::remove(grandparent->expressions.begin(),
                                                                   grandparent->expressions.end(),
                                                                   r->parent));
                        grandparent->expressions.push_back(query->whereClause);
                        query->parent->children.erase(std::remove(query->parent->children.begin(),
                                                                  query->parent->children.end(), query));
                        // Update preprocessor's state
                        this->qIndex.erase(query->id);
                    }
                }
            }
            else {
                leaderLogging("[ERROR] Unsupported EXISTS transformation.\n", true);
            }
        }
    }
    for (auto c : query->children)
        this->preprocess(c, db);
    if (query->sibling != nullptr)
        this->preprocess(query->sibling, db);
}



// Transforms the parsed query into a query graph
std::shared_ptr<QueryGraph> Preprocessor::transform(Parser::SQLQuery &query, Database &db) {
    // Collect query info
    auto qd = this->collectInfo(query);
    leaderLogging("[INFO] Collected " + std::to_string(this->qc) + " queries.\n");
    // Validate query and its sub-queries
    for (auto rit : this->qIndex) {
        this->resolveSchema(rit.second, db);
        this->validate(rit.second, db);
    }
    // Preprocess query
    this->preprocess(qd, db);
    qd->constructEdges(db);
    // Construct query graph
    qd->print();
    return qd;
}