#ifndef SECRECY_RULE_H
#define SECRECY_RULE_H

#include "plan.h"
#include "sql_expression.h"

/***
 * A plan transformation rule.
 */
class Rule {
protected:
    struct MatchInfo {
        std::shared_ptr<Plan> op1 = nullptr;
        std::shared_ptr<Plan> op2 = nullptr;
        // Additional info for Rules 1 and 2
        unsigned origin;
        // Additional info for Rule6
        std::vector<std::shared_ptr<SQLExpression>> predicates;
        std::shared_ptr<Expression> e1 = nullptr;
        std::shared_ptr<Expression> e2 = nullptr;
    };
    // Does pattern matching
    virtual std::shared_ptr<MatchInfo> match(std::shared_ptr<Plan> plan) = 0;
    // Applies the actual rule to create a new plan
    virtual std::shared_ptr<Plan> transform(std::shared_ptr<Plan> plan, std::shared_ptr<MatchInfo> info) = 0;
public:
    /***
     * Constructor
     */
    Rule() { };
    /***
     * Destructor
     */
    virtual ~Rule() { };

    /***
     * Applies the rule to the given plan.
     * @param plan - The input plan.
     * @return A new plan with lower cost or nullptr if the rule cannot be applied.
     */
    std::shared_ptr<Plan> apply(std::shared_ptr<Plan> plan);
};

// Applies semi-join transformation
class Rule1 : public Rule {
private:
    // Does pattern matching
    std::shared_ptr<Rule::MatchInfo> match(std::shared_ptr<Plan> plan);
    // Applies the actual rule to create a new plan
    std::shared_ptr<Plan> transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info);
public:
    /***
     * Constructor
     */
    Rule1() { };
    /***
     * Destructor
     */
    ~Rule1() { };
};

// Applies group-by decomposition
class Rule2 : public Rule {
private:
    // Does pattern matching
    std::shared_ptr<Rule::MatchInfo> match(std::shared_ptr<Plan> plan);
    // Applies the actual rule to create a new plan
    std::shared_ptr<Plan> transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info);
public:
    /***
     * Constructor
     */
    Rule2() { };
    /***
     * Destructor
     */
    ~Rule2() { };
};

// Applies select-distinct fusion
class Rule3 : public Rule {
private:
    // Does pattern matching
    std::shared_ptr<Rule::MatchInfo> match(std::shared_ptr<Plan> plan);
    // Applies the actual rule to create a new plan
    std::shared_ptr<Plan> transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info);
public:
    /***
     * Constructor
     */
    Rule3() { };
    /***
     * Destructor
     */
    ~Rule3() { };
};

// Applies self-join elimination
class Rule5 : public Rule {
private:
    // Does pattern matching
    std::shared_ptr<Rule::MatchInfo> match(std::shared_ptr<Plan> plan);
    // Applies the actual rule to create a new plan
    std::shared_ptr<Plan> transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info);
public:
    /***
     * Constructor
     */
    Rule5() { };
    /***
     * Destructor
     */
    ~Rule5() { };
};

// Applies distinct push-down
class Rule6 : public Rule {
private:
    // Does pattern matching
    std::shared_ptr<Rule::MatchInfo> match(std::shared_ptr<Plan> plan);
    // Applies the actual rule to create a new plan
    std::shared_ptr<Plan> transform(std::shared_ptr<Plan> plan, std::shared_ptr<Rule::MatchInfo> info);
public:
    /***
     * Constructor
     */
    Rule6() { };
    /***
     * Destructor
     */
    ~Rule6() { };
};

// Applies sort push-down
static std::shared_ptr<Plan> rule4(const std::shared_ptr<const Plan> plan) {
    return nullptr;
}

#endif //SECRECY_RULE_H
