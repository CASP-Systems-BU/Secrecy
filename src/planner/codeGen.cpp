#include "../../include/planner/codeGen.h"

// Generates physical operators for the given logical plan and database
std::shared_ptr<Operator> CodeGen::generatePhysicalPlan(std::vector<std::shared_ptr<Plan>> &plans, Database *db) {
    this->db = db;
    std::shared_ptr<Operator> root = nullptr;
    for (auto plan : plans) {
        std::vector<std::shared_ptr<Plan>> operators;
        plan->collectOperators(operators);

        std::vector<std::shared_ptr<Operator>> plan_ops;
        plan_ops.resize(operators.size(), nullptr);
        auto relations = db->getRelations();
        unsigned cnt = this->generateOperators(operators, relations, plan_ops);

        leaderLogging("[INFO] Generated " + std::to_string(cnt) + " physical operators.\n");

        root = plan_ops.back();
    }
    assert(root != nullptr);
    return root;
}

// Retrieves relation id
unsigned CodeGen::getRelationInd(const std::vector<std::shared_ptr<Relation>> &relations, std::string name) {
    for (unsigned i = 0; i < relations.size(); i++) {
        if (relations[i]->getName() == name) {
            return i;
        }
    }
    leaderLogging("[ERROR] Relation name '" + name + "' not found.\n", true);
    return relations.size();
}

// Retrieves index of input operator
unsigned CodeGen::getInputOp(std::shared_ptr<Plan> inputOp, const std::vector<std::shared_ptr<Plan>> &operators) {
    for (unsigned i = 0; i < operators.size(); i++) {
        if (operators[i] == inputOp) {
            return i;
        }
    }
    leaderLogging("[ERROR] Input operator not found.\n", true);
    return operators.size();
}

// Generates physical operators
unsigned CodeGen::generateOperators(std::vector<std::shared_ptr<Plan>> &operators,
                                    const std::vector<std::shared_ptr<Relation>> &relations,
                                    std::vector<std::shared_ptr<Operator>> &plan_ops) {
    unsigned op_i = 0;

    for (auto op: operators) {
        switch (op->op) {
            case Plan::Scan:
            { // Create scan operator
                std::shared_ptr<ScanExpr> s_expr = std::dynamic_pointer_cast<ScanExpr>(op->exp);
                auto rName = s_expr->rel->name;
                unsigned rel_i = getRelationInd(relations, rName);
                leaderLogging("[INFO] Creating Scan on relation " + rName + "...\n");
                plan_ops[op_i++] = std::make_shared<ScanOperator>(relations[rel_i],
                                                                  relations[rel_i]->getSize());
                break;
            }
            case Plan::Filter:
            { // Create filter operator
                unsigned ind = getInputOp(op->left, operators);
                auto input = plan_ops[ind];
                leaderLogging("[INFO] Creating Filter...\n");
                std::shared_ptr<FilterExpr> f_expr = std::dynamic_pointer_cast<FilterExpr>(op->exp);
                assert(f_expr != nullptr);
                plan_ops[op_i++] = std::make_shared<SelectOperator>(input, f_expr->selPredicate,
                                                                    input->batchSize);
                break;
            }
            case Plan::JoinFilter:
            { // Create filter operator
                unsigned ind = getInputOp(op->left, operators);
                auto input = plan_ops[ind];
                leaderLogging("[INFO] Creating JoinFilter...\n");
                auto j_expr = std::dynamic_pointer_cast<JoinExpr>(op->exp);
                assert(j_expr != nullptr);
                plan_ops[op_i++] = std::make_shared<SelectOperator>(input, j_expr->predicate,
                                                                    input->batchSize);
                break;
            }
            case Plan::NestedLoopJoin:
            { // Create join operator
                std::shared_ptr<JoinExpr> j_expr = std::dynamic_pointer_cast<JoinExpr>(op->exp);
                unsigned ind1 = getInputOp(op->left, operators);
                unsigned ind2 = getInputOp(op->right, operators);
                auto left = plan_ops[ind1];
                auto right = plan_ops[ind2];
                if (j_expr->jType == JoinExpr::Semi) { // Semi-join
                    leaderLogging("[INFO] Creating SemiJoin...\n");
                    if (j_expr->agg.size() == 0) {
                        plan_ops[op_i++] = std::make_shared<SemiJoinOperator>(left, right,
                                                                              j_expr->predicate,
                                                                              left->batchSize);
                    }
                    else { // Has partial aggregation
                        plan_ops[op_i++] = std::make_shared<SemiJoinOperator>(left, right,
                                                                              j_expr->predicate,
                                                                              j_expr->agg,
                                                                              left->batchSize);
                    }
                }
                else { // Inner or outer join
                    leaderLogging("[INFO] Creating Join...\n");
                    plan_ops[op_i++] = std::make_shared<JoinOperator>(left, right,
                                                                      j_expr->predicate,
                                                                      left->batchSize,
                                                                      right->batchSize);
                }
                break;
            }
            case Plan::GroupBy:
            { // Create group-by operator
                unsigned ind = getInputOp(op->left, operators);
                auto input = plan_ops[ind];
                leaderLogging("[INFO] Creating GroupBy...\n");
                std::shared_ptr<GroupByExpr> g_expr = std::dynamic_pointer_cast<GroupByExpr>(op->exp);
                int num_aggr = g_expr->aggs.size();  // TODO: handle more than one aggregations
                assert(!g_expr->aggs.empty());
                GroupByType aggr_type = GroupByNameType[g_expr->aggs[0]->name];
                std::vector<std::shared_ptr<Expression>> result_expressions;
//                result_expressions.push_back(std::make_shared<AttributeExpr>("R*S", "SEL"));
//                result_expressions.push_back(std::make_shared<AttributeExpr>("R*S", "[SEL]"));
                plan_ops[op_i++] = std::make_shared<GroupByOperator>(input,
                                                                   aggr_type, g_expr->keys,
                                                                   result_expressions,
                                                                   true,
                                                                   input->batchSize);


                break;
            }
            case Plan::OddEvenMerge:
            { // Create second phase of group-by operator
                unsigned ind = getInputOp(op->left, operators);
                auto input = plan_ops[ind];
                leaderLogging("[INFO] Creating OddEvenAgg...\n");
                std::shared_ptr<GroupByExpr> g_expr = std::dynamic_pointer_cast<GroupByExpr>(op->exp);
                GroupByType aggr_type = GroupByNameType[g_expr->aggs[0]->name];
                plan_ops[op_i++] = std::make_shared<GroupByOperator>(input,
                                                                     aggr_type, g_expr->keys,
                                                                     false, input->batchSize,
                                                                     true);
                break;
            }
            case Plan::OrderBy:
            { // Create sort operator
                unsigned ind = getInputOp(op->left, operators);
                auto input = plan_ops[ind];
                leaderLogging("[INFO] Creating OrderBy...\n");
                std::shared_ptr<OrderByExpr> o_expr = std::dynamic_pointer_cast<OrderByExpr>(op->exp);
                std::vector<bool> asc = o_expr->asc;
                plan_ops[op_i++] = std::make_shared<SortOperator>(input, o_expr->keys, asc,
                                                                  input->batchSize);
                break;
            }
            case Plan::Distinct:
            { // Create distinct operator
                unsigned ind = getInputOp(op->left, operators);
                auto input = plan_ops[ind];
                leaderLogging("[INFO] Creating Distinct...\n");
                std::shared_ptr<DistinctExpr> d_expr = std::dynamic_pointer_cast<DistinctExpr>(op->exp);
                plan_ops[op_i++] = std::make_shared<DistinctOperator>(input, d_expr->keys,
                                                                      true, input->batchSize);
                break;
            }
            case Plan::AdjacentEq:
            { // Create second phase of distinct operator
                unsigned ind = getInputOp(op->left, operators);
                auto input = plan_ops[ind];
                leaderLogging("[INFO] Creating AdjacentEq...\n");
                std::shared_ptr<DistinctExpr> d_expr = std::dynamic_pointer_cast<DistinctExpr>(op->exp);
                plan_ops[op_i++] = std::make_shared<DistinctOperator>(input, d_expr->keys,
                                                                      false, input->batchSize);
                break;
            }
            case Plan::TableFunction:
            {
                unsigned ind = getInputOp(op->left, operators);
                auto input = plan_ops[ind];
                leaderLogging("[INFO] Creating TableFunction...\n");
                std::shared_ptr<TableFunctionExpr> f_expr = std::dynamic_pointer_cast<TableFunctionExpr>(op->exp);
                assert(f_expr);
                FunctionType func_type = FunctionNameType[f_expr->exp->name];
                plan_ops[op_i++] = std::make_shared<TableFunctionOperator>(input, func_type, f_expr->exp);
//                leaderLogging("[ERROR] Table function not supported yet.\n", true);
                break;
            }
            case Plan::Mask:
            { // Create mask operator
                unsigned ind = getInputOp(op->left, operators);
                auto input = plan_ops[ind];
                leaderLogging("[INFO] Creating Mask...\n");
                plan_ops[op_i++] = std::make_shared<MaskOperator>(input, input->batchSize);
                break;
            }
            default:
            {
                std::cout << op->exp->toString() << std::endl;
                leaderLogging("[ERROR] Unknown operator.\n", true);
                break;
            }
        }
    }
    return op_i;
}