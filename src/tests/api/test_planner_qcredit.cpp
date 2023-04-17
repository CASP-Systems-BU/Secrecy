#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"
#include <map>
#include <set>

// TODO: Increase the number of rows.
bool checkCorrectness(Data **r_data, RowData res_open,
                      int uid_ind, int year_ind, int score_ind,
                      int sel_ind,
                      int size) {

    if (get_rank() == 0) {

        std::map<Data, Data> credit_min;
        std::map<Data, Data> credit_max;
        for (int i = 0; i < size; ++i) {
            if(r_data[i][year_ind] == 50) {
                if (credit_min.find(r_data[i][uid_ind]) != credit_min.end()) {
                    credit_min[r_data[i][uid_ind]] = std::min(r_data[i][score_ind], credit_min[r_data[i][uid_ind]]);
                    credit_max[r_data[i][uid_ind]] = std::max(r_data[i][score_ind], credit_max[r_data[i][uid_ind]]);
                } else {
                    credit_min[r_data[i][uid_ind]] = r_data[i][score_ind];
                    credit_max[r_data[i][uid_ind]] = r_data[i][score_ind];
                }
            }
        }

        std::set<Data> real_set;
        for(auto it = credit_max.begin(); it != credit_max.end(); ++it){
            if ((it->second - credit_min[it->first]) > 10){
                real_set.insert(it->first);
            }
        }

        auto res_open_ptr = res_open.get();
        std::set<Data> res_set;
        for (int i = 0; i < size; ++i) {
            if (res_open_ptr[sel_ind][i]) {
                res_set.insert(res_open_ptr[uid_ind][i]);
            }
        }

        printf("Users count %ld==%ld\n", real_set.size(), res_set.size());
        return real_set == res_set;
    } else {
        return true;
    }
}

void planner_qcredit_test() {
    // Table schema
    int credits_size = 256;
    std::vector<std::string> credits_schema = {"[uid]", "[year]",
                                               "[year-2019]", "[2019-year]",
                                               "[score]", "[score+threshold]",
                                               "[SEL]", "SEL"};
    Data **credits_data = getRandomData(credits_size, credits_schema.size());
    for (int i = 0; i < credits_size; ++i) {
        credits_data[i][2] = credits_data[i][1] - 50;
        credits_data[i][3] = 50 - credits_data[i][1];
        credits_data[i][5] = credits_data[i][1] + (credits_data[i][4] % 10);
    }
    Relation credits = getTestRelation("CREDITS", credits_size, credits_schema, credits_data);

    //exchange seeds
    exchange_rsz_seeds();

    auto uid_expr = Expression("CREDITS", "[uid]");
    auto score_expr = Expression("CREDITS", "[score]");
    auto threshold_expr = Expression("CREDITS", "[score+threshold]");
    auto sel_expr = Expression("CREDITS", "[SEL]");


    // Query
    // Selection
    credits["[SEL]"] = !((credits["[year-2019]"] < 0) ^ (credits["[2019-year]"] < 0))
                        & (credits[score_expr] > credits[threshold_expr]);
    credits.bSelectionUsed = true;

    ScanOperator scan(credits, 4);
    GroupByOperator groupBy(scan, GroupByType::MIN_MAX, {sel_expr, uid_expr}, {score_expr, threshold_expr}, true, 4);
    SortOperator sort(groupBy, {uid_expr}, {1}, 4);
    sort.execute();

    auto res_open = sort.output.openRelation();
    bool test_res = checkCorrectness(credits_data, res_open,
                                     0, 1, 4,6,
                                     credits_size);

    assert(test_res == true);

    // free temp data
    free(credits_data);
    credits.freeRelation();
}


int main(int argc, char **argv) {

    // initialize communication
    init(argc, argv);
    init_sharing();

    // Testing and timing the query
    timeTest(planner_qcredit_test, "planner_qcredit_test");

    // tear down communication
    MPI_Finalize();
    return 0;
}