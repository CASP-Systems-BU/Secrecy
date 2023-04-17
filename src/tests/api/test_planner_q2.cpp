#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"
#include <map>

bool checkQ2Correctness(Data **r_data, Data **s_data,
                        RowData res_open,
                        int r_id, int r_ak, int s_id,
                        int res_ak, int res_sel, int res_sel_a,
                        int rows_1, int rows_2, int rows_3) {

    if (get_rank() == 0) {
        std::map<Data, long long> real_count;
        for (int i_1 = 0; i_1 < rows_1; ++i_1) {
            for (int i_2 = 0; i_2 < rows_2; ++i_2) {
                if (r_data[i_1][r_id] == s_data[i_2][s_id]) {
                    real_count[r_data[i_1][r_ak]]++;
                    // printf("%d,%d:\t%lld==%lld\tcount(%lld)=%lld\n", i_1, i_2, r_data[i_1][r_id], s_data[i_2][s_id], r_data[i_1][r_ak], real_count[r_data[i_1][r_ak]]);
                }
            }
        }

        auto res_open_ptr = res_open.get();
        std::map<Data, long long> res_count;
        for (int i = 0; i < rows_3; ++i) {
            if (res_open_ptr[res_ak][i] != -1) {
                res_count[res_open_ptr[res_ak][i]] = res_open_ptr[res_sel_a][i];
                // printf("%d:\t\tcount(%lld)=%lld\n", i, res_open_ptr[res_ak][i], res_open_ptr[res_sel_a][i]);
            }
        }

        // printf("Groups count %ld==%ld\n", res_count.size(), real_count.size());
        // for(auto it = res_count.cbegin(); it != res_count.cend(); ++it){
        //     std::cout << it->first << "\t" << it->second << "==" << real_count[it->first] << "\n";
        // }

        printf("Groups count %ld==%ld\n", res_count.size(), real_count.size());

        return real_count == res_count;
    } else {
        return true;
    }
}

void planner_q2_test() {
    // Table schema
    int r_rows = 128;
    std::vector<std::string> r_schema = {"[id]", "[ak]", "[SEL]"};
    Data **r_data = getRandomData(r_rows, r_schema.size());
    // adding selection ones in arth and b
    for (int i = 0; i < r_rows; ++i) {
        r_data[i][2] = 1;
    }
    Relation r = getTestRelation("R", r_rows, r_schema, r_data);
    r.aSelectionUsed = false;
    r.bSelectionUsed = true;


    // S table
    int s_rows = 128;
    std::vector<std::string> s_schema = {"[id]", "[SEL]"};
    Data **s_data = getRandomData(s_rows, s_schema.size());
    Relation s = getTestRelation("S", s_rows, s_schema, s_data);

    //exchange seeds
    exchange_rsz_seeds();

    // SELECT R.ak , COUNT(*)
    // FROM R, S
    // WHERE R.id = S.id
    // GROUP BY R.ak

    // Query
    ScanOperator r_Scan(r);
    ScanOperator s_Scan(s);
    Expression join_expr(ExpressionType::Equal, Expression("R", "[id]"), Expression("S", "[id]"));
    JoinOperator rs_Join(r_Scan, s_Scan, join_expr);
    Expression star(ExpressionType::Star);
    GroupByOperator rs_group_by(rs_Join, GroupByType::COUNT, {Expression("R*S", "[SEL]"), Expression("R", "[ak]")}, {Expression("R*S", "SEL"), Expression("R*S", "[SEL]")}, {GroupByType::COUNT}, {star});

    // Execution
    rs_group_by.execute();

    auto res_open = rs_group_by.output.openRelation();

    const int r_id = r["[id]"].getIndex();
    const int r_ak = r["[ak]"].getIndex();
    const int s_id = s["[id]"].getIndex();
    const int res_ak = rs_group_by.output["R.[ak]"].getIndex();
    const int res_sel = rs_group_by.output["[SEL]"].getIndex();
    const int res_sel_a = rs_group_by.output["SEL"].getIndex();
    const int rows_3 = rs_group_by.output.getSize();

    bool test_res = checkQ2Correctness(r_data, s_data,
                                       res_open,
                                       r_id, r_ak, s_id,
                                       res_ak, res_sel, res_sel_a,
                                       r_rows, s_rows, rows_3);

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
    timeTest(planner_q2_test, "planner_q2_test");

    // tear down communication
    MPI_Finalize();
    return 0;
}