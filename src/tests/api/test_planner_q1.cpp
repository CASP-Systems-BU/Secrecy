#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"
#include <set>

bool checkQ1Correctness(Data **r_data, Data **s_data,
                        RowData res_open,
                        int r_id, int s_id, int res_id, int res_sel,
                        int rows_1, int rows_2, int rows_3) {

    if (get_rank() == 0) {
        // First table set
        std::set<Data> r_b_data_set;
        for (int i = 0; i < rows_1; ++i) {
            r_b_data_set.insert(r_data[i][r_id]);
        }

        std::set<Data> s_b_data_set;
        for (int i = 0; i < rows_2; ++i) {
            s_b_data_set.insert(s_data[i][s_id]);
        }

        std::set<Data> intersection;
        set_intersection(r_b_data_set.begin(), r_b_data_set.end(), s_b_data_set.begin(), s_b_data_set.end(),
                         std::inserter(intersection, intersection.begin()));

        // res table set
        auto res_open_ptr = res_open.get();
        std::set<Data> res_set;
        for (int i = 0; i < rows_3; ++i) {
            if (res_open_ptr[res_sel][i]) {
                res_set.insert(res_open_ptr[res_id][i]);
            }
        }
        printf("Distinct count %ld==%ld\n", intersection.size(), res_set.size());

        return intersection == res_set;
    } else {
        return true;
    }
}

void planner_q1_test() {
    // Table schema
    int r_rows = 128;
    std::vector<std::string> r_schema = {"[id]", "[SEL]"};
    Data **r_data = getRandomData(r_rows, r_schema.size());
    Relation r = getTestRelation("R", r_rows, r_schema, r_data);

    // S table
    int s_row = 128;
    std::vector<std::string> s_schema = {"[id]", "[SEL]"};
    Data **s_data = getRandomData(s_row, s_schema.size());
    Relation s = getTestRelation("S", s_row, s_schema, s_data);

    //exchange seeds
    exchange_rsz_seeds();

    // Query SQL:
    // SELECT DISTINCT R.id
    // FROM R, S
    // WHERE R.id = S.id

    // Query Baseline
    ScanOperator r_Scan(r, 4);
    ScanOperator s_Scan(s, 4);
    auto join_expr = std::make_shared<Expression>(Expression(Expression::Equal));
    join_expr->expressions.push_back(std::make_shared<AttributeExpr>(AttributeExpr("[id]", "R")));
    join_expr->expressions.push_back(std::make_shared<AttributeExpr>(AttributeExpr("[id]", "S")));
    JoinOperator rs_Join(r_Scan, s_Scan, join_expr, 4, 4);
    DistinctOperator rs_distinct(rs_Join, join_expr->expressions, true, 4);

    // Execution
    rs_distinct.execute();

    auto res_open = rs_distinct.output.openRelation();

    const int r_id = r["[id]"].getIndex();
    const int s_id = s["[id]"].getIndex();
    const int res_id = rs_distinct.output["R.[id]"].getIndex();
    const int res_sel = rs_distinct.output["[SEL]"].getIndex();
    const int rows_3 = rs_distinct.output.getSize();

    bool test_res = checkQ1Correctness(r_data, s_data,
                                       res_open,
                                       r_id, s_id, res_id, res_sel,
                                       r_rows, s_row, rows_3);

    assert(test_res == true);

    // free temp data
    free(r_data);
    free(s_data);
    r.freeRelation();
    s.freeRelation();
}


int main(int argc, char **argv) {

    // initialize communication
    init(argc, argv);
    init_sharing();

    // Testing and timing the query
    timeTest(planner_q1_test, "planner_q1_test");

    // tear down communication
    MPI_Finalize();
    return 0;
}