#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"
#include <set>

bool checkQ3Correctness(Data **r_data, Data **s_data,
                        RowData res_open,
                        int r_ak, int s_ab, int r_id,
                        int res_id, int res_sel,
                        int rows_1, int rows_2, int rows_3) {

    if (get_rank() == 0) {
        // First table set
        std::set<Data> real_set;
        for (int i_1 = 0; i_1 < rows_1; ++i_1) {
            for (int i_2 = 0; i_2 < rows_2; ++i_2) {
                if (r_data[i_1][r_ak] == s_data[i_2][s_ab]) {
                    real_set.insert(r_data[i_1][r_id]);
                }
            }
        }

        // res table set
        auto res_open_ptr = res_open.get();
        std::set<Data> res_set;
        for (int i = 0; i < rows_3; ++i) {
            if (res_open_ptr[res_sel][i] != 0) {
                res_set.insert(res_open_ptr[res_id][i]);
            }
        }

        printf("Distinct count %ld==%ld\n", real_set.size(), res_set.size());

        return real_set == res_set;
    } else {
        return true;
    }
}

void planner_q3_test() {
    // Table schema
    int r_rows = 128;
    std::vector<std::string> r_schema = {"[id]", "[ak]", "[ab]", "SEL", "[SEL]"};
    Data **r_data = getRandomData(r_rows, r_schema.size());
    Relation r_1 = getTestRelation("R1", r_rows, r_schema, r_data);
    Relation r_2 = getTestRelation("R2", r_rows, r_schema, r_data);

    //exchange seeds
    exchange_rsz_seeds();

    // Query:
    // SELECT DISTINCT R.id
    // FROM R
    // WHERE R.ak = R.ab

    // Query Baseline
    ScanOperator r1_Scan(r_1);
    ScanOperator r2_Scan(r_2);
    Expression join_expr(ExpressionType::Equal, Expression("R1", "[ak]"), Expression("R2", "[ab]"));
    JoinOperator rs_Join(r1_Scan, r2_Scan, join_expr, 8, 8);
    DistinctOperator rs_distinct(rs_Join,
                                 {Expression("R1*R2", "[SEL]"), Expression("R1", "[id]"), Expression("R2", "[id]")}, 8);

    // Execution
    rs_distinct.execute();

    auto res_open = rs_distinct.output.openRelation();

    const int r_ak = r_1["[ak]"].getIndex();
    const int r_id = r_1["[id]"].getIndex();
    const int s_ab = r_1["[ab]"].getIndex();
    const int res_id = rs_distinct.output["R1.[id]"].getIndex();
    const int res_sel = rs_distinct.output["[SEL]"].getIndex();
    const int rows_3 = rs_distinct.output.getSize();

    bool test_res = checkQ3Correctness(r_data, r_data,
                                       res_open,
                                       r_ak, s_ab, r_id, res_id, res_sel,
                                       r_rows, r_rows, rows_3);

    assert(test_res == true);

    // free temp data
    free(r_data);
    r_1.freeRelation();
    r_2.freeRelation();
}


int main(int argc, char **argv) {

    // initialize communication
    init(argc, argv);
    init_sharing();

    // Testing and timing the query
    timeTest(planner_q3_test, "planner_q3_test");

    // tear down communication
    MPI_Finalize();
    return 0;
}