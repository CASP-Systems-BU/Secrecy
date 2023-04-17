#ifndef SECRECY_EXPRESSION_H
#define SECRECY_EXPRESSION_H

#include <string>
#include <vector>
#include <map>
#include <set>
#include <iostream>

/***
 * A Secrecy expression.
 */
class Expression : public std::enable_shared_from_this<Expression> {
public:
    /***
     * Expression types
     */
    enum ExpressionType {
        RelationRef,
        AttributeRef,
        DerivedAttributeRef,
        FunctionRef,
        GreaterThan,
        LessThan,
        Equal,
        GreaterThanEqual,
        LessThanEqual,
        NotEqual,
        Addition,
        Subtraction,
        Multiplication,
        Division,
        And,
        Or,
        Not,
        Star,
        In,
        Exists,
        LiteralInt,
        LiteralFloat,
        LiteralString,
        Unknown
    };
private:
    /***
     * Expression type-to-string mapping
     */
    const std::map<ExpressionType, std::string> ExpressionTypeName = {
            {ExpressionType::RelationRef,      "Relation"},
            {ExpressionType::AttributeRef,     "Attribute"},
            {ExpressionType::FunctionRef,      "Function"},
            {ExpressionType::GreaterThan,      "GreaterThan"},
            {ExpressionType::LessThan,         "LessThan"},
            {ExpressionType::Equal,            "Equal"},
            {ExpressionType::GreaterThanEqual, "GreaterThanEqual"},
            {ExpressionType::LessThanEqual,    "LessThanEqual"},
            {ExpressionType::NotEqual,         "NotEqual"},
            {ExpressionType::Addition,         "Addition"},
            {ExpressionType::Subtraction,      "Subtraction"},
            {ExpressionType::Multiplication,   "Multiplication"},
            {ExpressionType::Division,         "Division"},
            {ExpressionType::And,              "And"},
            {ExpressionType::Or,               "Or"},
            {ExpressionType::Not,              "Not"},
            {ExpressionType::Star,             "Star"},
            {ExpressionType::In,               "In"},
            {ExpressionType::Exists,           "Exists"},
            {ExpressionType::LiteralInt,       "LiteralInt"},
            {ExpressionType::LiteralFloat,     "LiteralFloat"},
            {ExpressionType::LiteralString,    "LiteralString"},
            {ExpressionType::Unknown,          "Unknown"}
    };
public:
    // Denotes whether the expression is an arithmetic share
    bool arithmeticShare;
    // Denotes whether `this` is a special expression that compares adjacent table rows
    bool adjRows = false;
    /***
     * Expression name (e.g., attribute's name, function's name, etc.)
     */
    std::string name;
    /***
     * The name or alias of the relation an attribute belongs to (for attribute expressions)
     */
    std::string relation;
    /***
     * The relation id (for relation expressions)
     */
    int id;
    /***
     * The id of the logical operator `this` expression corresponds to.
     */
    int opId;
    /***
     * A flag that denotes whether `this` expression is a relation that represents the output of a sub-query
     */
    bool subQuery = false;
    /***
     * A flag that denotes whether `this` expression is a relation that corresponds to a sub-query with a custom schema
     */
    bool customSchema = false;
    /***
    * The relation schema
    */
    std::vector<std::shared_ptr<Expression>> schema;
    /***
     * The relation cardinality (for relation expressions)
     */
    double cardinality = 0;
    /***
     * A flag that denotes whether `this` expression is a function applied to distinct values, e.g., COUNT(DISTINCT id)
     */
    bool distinct = false;
    /***
     * String value (for string constants)
     */
    std::string sval;
    /***
     * Integer value (for integer constants)
     */
    int64_t ival;
    /***
     * Double value (for floating point constants)
     */
    double fval;
    /***
     * The id of the query this expression belongs to
     */
    int qID = -1;
    /***
     * Expression alias
     */
    std::string alias;
    /***
     * Expression type
     */
    ExpressionType expressionType;
    /***
     * Subexpressions (can be more than two, e.g., in ternary operations)
     */
    std::vector<std::shared_ptr<Expression>> expressions;
    /***
     * Parent expression
     */
    std::shared_ptr<Expression> parent = nullptr;

    /***
     * Constructor
     * @param eType - The type of the expression.
     * @param alias - The expression alias (if any).
     */
    Expression(const ExpressionType &eType, const char* alias=nullptr, std::shared_ptr<Expression> parent=nullptr) {
        this->expressionType = eType;
        this->parent = parent;
        if (alias != nullptr) {
            this->alias = std::string(alias);
        }
    }

    /***
     * Checks if two expressions are the same.
     * @param other - The expression to compare with `this` one.
     * @return True if `this` and `other` correspond to the same expression, False otherwise.
     */
    bool sameAs(std::shared_ptr<Expression> other) {
        if (expressionType == other->expressionType) {
            if ((name == other->name) && (relation == other->relation)) {
                return true;
            }
        }
        return false;
    }

    /***
     * Checks if a sub-query is contained in `this` expression.
     * @param sub_query_name - The name of the relation that corresponds to the sub-query result.
     * @return True if the expression contains the sub-query with name `sub_query_name`, False otherwise.
     */
    bool contains(std::string &sub_query_name) const {
        if (this->hasType(ExpressionType::RelationRef)) {
            if (this->toString() == sub_query_name)
                return true;
        }
        // Check sub-expressions
        for (auto e : expressions) {
            if (e->contains(sub_query_name))
                return true;
        }
        return false;
    }

    /***
     * Adds a sub-expression
     * @param expression - The sub-expression to add.
     */
    void addExpression(std::shared_ptr<Expression> expression) {
        this->expressions.push_back(expression);
    }

    /***
     * Checks if the expression is of type `expressionType`.
     * @param expressionType - The expression type to check.
     * @return True if the expression is of type `expressionType`, False otherwise.
     */
    bool hasType(const ExpressionType expressionType) const {
        return this->expressionType == expressionType;
    }

    /***
     * @return True if `this` is a logical expression, False otherwise.
     */
    bool isLogical() const {
        bool res = expressionType == And || 
                    expressionType == Or || 
                    expressionType == Not;
        return res;
    }

    /***
     * @return True if `this` is an arithmetic expression, False otherwise.
     */
    bool isArithmetic() const {
        bool res = expressionType == Addition || 
                    expressionType == Subtraction || 
                    expressionType == Multiplication ||
                    expressionType == Division;
        return res;
    }

    /***
     * @return True if `this` expression is a comparison, False otherwise.
     */
    bool isComparison() const {
        bool res = expressionType == GreaterThan ||
                   expressionType == LessThan ||
                   expressionType == Equal ||
                   expressionType == GreaterThanEqual ||
                   expressionType == LessThanEqual ||
                   expressionType == NotEqual;
        return res;
    }

    /***
     * Traverses the logical expression tree and collects individual predicates.
     * @param predicates - The vector of predicates to populate.
     */
    void flatten(std::vector<std::shared_ptr<Expression>> &predicates) {
        if (this->isLogical()) {
            for (auto exp: expressions)
                exp->flatten(predicates);
        }
        else {
            predicates.push_back(shared_from_this());
        }
    }

    /***
     * @return True if `this` expression has the form 'Atom_1 AND Atom_2 AND ... AND Atom_k', False otherwise.
     */
    bool isConjunctive() const {
        bool r = true;
        if (this->isLogical()) {
            if (expressionType != And)
                return false;
            for (auto e: expressions)
                r &= e->isConjunctive();
        }
        return r;
    }

    /***
     * @return The root of the tree `this` expression belongs to.
     */
    std::shared_ptr<Expression> findRoot() {
        if (parent)
            return parent->findRoot();
        return shared_from_this();
    }

    /***
     * Collects attributes in `this` expression.
     * @param attributes - The vector of attributes to populate.
     */
    void collectAttributes(std::vector<std::shared_ptr<Expression>> &attributes) {
        if (expressionType != ExpressionType::AttributeRef) {
            for (auto exp: expressions)
                exp->collectAttributes(attributes);
        }
        else {
            attributes.push_back(shared_from_this());
        }
    }

    /***
     * Constructs and returns a string representation of `this` expression.
     * @return A string representation of `this` expression.
     */
    virtual std::string toString() const {
        std::string str;
        switch (this->expressionType) {
            case ExpressionType::Star:
                str = "*";
                break;
            default:  // Continue recursively
                str = ExpressionTypeName.at(expressionType) + "[";
                for (auto expr: expressions) {
                    str += expr->toString();
                    str += ", ";
                }
                if (expressions.size() > 0) { // Remove last comma
                    str = str.substr(0, str.size() - 2);
                }
                str += "]";
        }
        return str;
    }

    /***
     * Retrieves the temporary relation expression given a name or alias.
     * @param rel_id - The name of alias of the temporary relation.
     * @return The temporary relation expression.
     */
    std::shared_ptr<Expression> getRelation(const std::string &rel_id) {
        if (this->hasType(ExpressionType::RelationRef)) {
            if (this->toString() == rel_id)
                return shared_from_this();
        }
        // Check sub-expressions
        for (auto e : expressions) {
            if (auto r = e->getRelation(rel_id))
                return r;
        }
        return nullptr;
    }

    /***
     * Retrieves all temporary relation expressions that appear in the query.
     * @param relations - The vector of relations to populate.
     */
    void collectTempRelations(std::vector<std::shared_ptr<Expression>> &relations) {
        if (expressionType != ExpressionType::RelationRef) {
            for (auto exp: expressions)
                exp->collectTempRelations(relations);
        }
        else {
            relations.push_back(shared_from_this());
        }
    }

    /***
     * Set the sharing types for `this` expression and all its subexpressions.
     * @param sTypes - A mapping from attribute names to boolean (True) or Arithmetic (False) share types
     */
    void setShareTypes(std::map<std::string, bool> &sTypes) {
        if (expressionType==AttributeRef) {
            auto it = sTypes.find(name);
            if (it != sTypes.end()) { // If it's a base attribute
                if (it->second) {
                    name = "[" + name + "]";
                    arithmeticShare = false;
                }
                else {
                    arithmeticShare = true;
                }
            }
        }
        else {
            for (auto e : expressions)
                e->setShareTypes(sTypes);
        }
    }
};

/***
 * A named expression.
 */
class NamedExpr : public Expression {
public:
    /***
     * Constructor.
     * @param eType - The type of the named expression.
     * @param name - The expression name.
     * @param alias - The expression alias (if any).
     */
    NamedExpr(ExpressionType eType,
              const char* name,
              const char* alias,
              std::shared_ptr<Expression> parent=nullptr) : Expression(eType, alias, parent)
    {
        assert(name != nullptr);
        this->name = std::string(name);
    }

    /***
     * @return The expression name.
     */
    virtual std::string getIdentifier() const {
        return name+"#"+alias;
    }

    /***
     * @return The expression name.
     */
    virtual std::string toString() const override {
        if (this->alias.empty())
            return name;
        else
            return name + alias;
    }
};

/***
 * A relation attribute.
 */
class AttributeExpr : public NamedExpr {
public:
    /***
     * Constructor.
     * @param name - The attribute name
     * @param relName - The name of the relation the attribute belongs to.
     * @param alias - The attribute alias (if any).
     */
    AttributeExpr(const char* name,
                  const char* relId,
                  const char* alias=nullptr,
                  std::shared_ptr<Expression> parent=nullptr) : NamedExpr(ExpressionType::AttributeRef, name,
                                                                          alias, parent)
    {
        if (relId != nullptr) {
            this->relation = std::string(relId);
        }
    }

    /***
     * Shallow copy constructor
     * @param other - The function expression to copy from.
     */
    AttributeExpr(std::shared_ptr<const AttributeExpr> other) : NamedExpr(ExpressionType::AttributeRef,
                                                                          "", "", nullptr)
    {
        this->name = other->name;
        this->alias = other->alias;
        this->expressions = other->expressions;
        this->opId = other->opId;
    }

    /***
     * Constructs and returns an attribute identifier of the form 'relation_name.attribute_name'.
     * @return The attribute identifier.
     */
    std::string getIdentifier() const override {
        return relation + "." + name;
    }

    /***
     * Constructs and returns a string representation of `this` attribute.
     * @return The string representation of `this` attribute.
     */
    std::string toString() const override {
        if (this->alias.empty())
            return getIdentifier();
        else
            return getIdentifier() + " as " + alias;
    }
};

/***
 * A base relation or a sub-query output.
 */
class RelationExpr : public NamedExpr {
public:
    /***
     * Constructor.
     * @param name - The relation name.
     * @param alias - The relation alias (if any).
     * @param schema - The custom relation schema (if any).
     * @param subQuery - A flag that denotes whether the relation is the output of a sub-query (default: False).
     */
    RelationExpr(const char* name,
                 const char* alias,
                 bool subQuery=false,
                 const std::vector<char*>* schema=nullptr,
                 std::shared_ptr<Expression> parent=nullptr) : NamedExpr(ExpressionType::RelationRef, name,
                                                                         alias, parent)
    {
        this->subQuery = subQuery;
        if (schema != nullptr) {
            this->customSchema = true;
            for (auto c : *schema) {
                auto n = alias != nullptr ? alias : name;
                auto att = std::make_shared<AttributeExpr>(AttributeExpr(c, n));
                this->schema.push_back(att);
            }
        }
    }

    /***
     * @return The relation name.
     */
    std::string toString() const override {
        if (this->alias.empty())
            return name;
        else
            return name + " as " + alias;
    }

    bool identifiedBy(const std::string &name_or_alias) {
        if (this->alias.empty() && this->name==name_or_alias)
            return true;
        else if (!this->alias.empty() && this->alias==name_or_alias)
            return true;
        return false;
    }

    /***
     * Checks if `this` relation contains the given attribute.
     * @param attribute - The attribute to search for.
     * @return True if the attribute exists in `this` relation, False otherwise.
     */
    bool contains(const std::shared_ptr<AttributeExpr> attribute) {
        if (!alias.empty() && (attribute->relation != alias))
            return false;
        if (alias.empty() && attribute->relation != name)
            return false;
        for (auto c : schema)
            if (c->name == attribute->name)
                return true;
        return false;
    }
};

/***
 * A function, e.g., COUNT, SUM, MIN, MAX, AVG, etc.
 */
class FunctionExpr : public NamedExpr {
public:
    /***
     * Constructor.
     * @param name
     * @param alias
     * @param distinct
     */
    FunctionExpr(const char* name,
                 const char* alias,
                 bool distinct=false,
                 std::shared_ptr<Expression> parent=nullptr) : NamedExpr(ExpressionType::FunctionRef, name,
                                                                         alias, parent)
    {
        this->distinct = distinct;
        // Convert function name to uppercase
        std::transform(this->name.begin(), this->name.end(),this->name.begin(), ::toupper);
    }

    /***
     * Shallow copy constructor
     * @param other - The function expression to copy from.
     */
    FunctionExpr(std::shared_ptr<const FunctionExpr> other) : NamedExpr(ExpressionType::FunctionRef,
                                                                        "", "", other->parent)
    {
        this->name = other->name;
        this->alias = other->alias;
        this->expressions = other->expressions;
        this->distinct = other->distinct;
        this->opId = other->opId;
    }

    /***
     * Constructs and returns a string representation of `this` function application.
     * @return The string representation of `this` function application.
     */
    std::string toString() const override {
        assert(expressions.size()>0);
        std::string str = name + "(";
        str += this->distinct ? "DISTINCT " : "";
        auto e = expressions[0];
        str += e->toString() + ")";
        str += this->alias.empty() ? "" : " as " + this->alias;
        return str;
    }
};

/***
 * A constant, i.e., a string literal, an integer, or a floating point number.
 */
class ConstantExpr : public Expression {
public:
    /***
     * String literal constructor.
     * @param sval
     * @param alias
     */
    ConstantExpr(const char* sval, const char* alias=nullptr,
                 std::shared_ptr<Expression> parent=nullptr) : Expression(ExpressionType::LiteralString, alias,
                                                                          parent)
    {
        assert(sval != nullptr);
        this->sval = std::string(sval);
    }

    /***
     * Integer literal constructor.
     * @param ival
     * @param alias
     */
    ConstantExpr(int64_t ival, const char* alias=nullptr,
                 std::shared_ptr<Expression> parent=nullptr) : Expression(ExpressionType::LiteralInt, alias,
                                                                          parent)
    {
        this->ival = ival;
    }

    /***
     * Float literal constructor.
     * @param fval
     * @param alias
     */
    ConstantExpr(double fval, const char* alias=nullptr,
                 std::shared_ptr<Expression> parent=nullptr) : Expression(ExpressionType::LiteralFloat, alias,
                                                                          parent)
    {
        this->fval = fval;
    }

    /***
     * Constructs and returns a string representation of `this` constant.
     * @return The string representation of `this` constant.
     */
    std::string toString() const override {
        switch (expressionType) {
            case LiteralString:
                return "'"+sval+"'";
            case LiteralInt:
                return std::to_string(ival);
            case LiteralFloat:
                return std::to_string(fval);
            default:
                std::cout << "[ERROR] Unsupported constant.\n" << std::endl;
                exit(-1);
        }
    }
};

/***
 * An attribute generated during query processing.
 */
class DerivedAttributeExpr : public NamedExpr {
public:
    /***
     * Constructor.
     * @param name - The name of the derived attribute.
     */
    DerivedAttributeExpr(const char *name, int opId=-1) : NamedExpr(ExpressionType::DerivedAttributeRef, name,
                                                                     nullptr, nullptr)
    {
        this->opId = opId;
    }

    /***
     * Constructs and returns a string representation of `this` attribute.
     * @return The string representation of `this` attribute.
     */
    std::string toString() const override {
        return relation.empty() ? name : relation + "." + name;
    }
};

#endif //SECRECY_EXPRESSION_H