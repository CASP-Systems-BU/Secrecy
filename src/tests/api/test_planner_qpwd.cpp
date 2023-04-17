#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"
#include <map>

bool checkCorrectness(Data **r_data, RowData res_open,
                      int pid_ind, int password_ind,
                      int sel_ind, int count_ind,
                      int size) {

    if (get_rank() == 0) {

        std::map<std::pair<Data, Data>, Data> all_count;
        std::map<std::pair<Data, Data>, Data> chosen_count;

        for (int i = 0; i < size; ++i) {
            all_count[{r_data[i][pid_ind], r_data[i][password_ind]}]++;
        }

        for (auto it = all_count.begin(); it != all_count.end(); ++it) {
            if (it->second > 1) {
                chosen_count[it->first] = it->second;
            }
        }

        auto res_open_ptr = res_open.get();
        std::map<std::pair<Data, Data>, Data> res_count;
        

        for (int i = 0; i < size; ++i) {
            if (res_open_ptr[sel_ind][i] && res_open_ptr[pid_ind][i] > 0) {
                res_count[{res_open_ptr[pid_ind][i], res_open_ptr[password_ind][i]}] = res_open_ptr[count_ind][i];
            }
        }

        printf("Duplicates count %ld==%ld\n", chosen_count.size(), res_count.size());

        return chosen_count == res_count;
    } else {
        return true;
    }
}

void planner_qpwd_test() {
    // Table schema
    int passwords_size = 1024;
    std::vector<std::string> passwords_schema = {"[uid]", "[pwd]", "[count]", "[SEL]", "SEL", "count"};
    Data **passwords_data = getRandomData(passwords_size, passwords_schema.size());
    for (int i = 0; i < passwords_size; ++i) {
        passwords_data[i][3] = 1;
        passwords_data[i][4] = 1;
    }
    Relation passwords = getTestRelation("PASSWORDS", passwords_size, passwords_schema, passwords_data);
    passwords.bSelectionUsed = true;

    //exchange seeds
    exchange_rsz_seeds();

    auto uid_expr = Expression("PASSWORDS", "[uid]");
    auto pwd_expr = Expression("PASSWORDS", "[pwd]");
    auto sel_b = Expression("PASSWORDS", "[SEL]");
    auto sel_a = Expression("PASSWORDS", "SEL");
    auto count_b = Expression("PASSWORDS", "[count]");
    auto count_a = Expression("PASSWORDS", "count");

    // Query
    ScanOperator scan(passwords, 4);
    SortOperator sort(scan, {uid_expr, pwd_expr}, {1, 1}, 4);
    Expression star(ExpressionType::Star);
    GroupByOperator groupBy(sort, GroupByType::COUNT, {uid_expr, pwd_expr}, {count_a, sel_b}, {GroupByType::COUNT}, {star}, false, 4);
    ConversionOperator conversion(groupBy, {count_a}, {count_b});
    conversion.execute();

    Relation passwords_temp = conversion.output;
    passwords_temp[sel_b] = passwords_temp[sel_b] & (passwords_temp[count_b] > 1);

    // TODO: create mask
    ScanOperator scan2(passwords_temp, 4);
    SortOperator sort_2(scan2, {uid_expr}, {1}, 4);
    sort_2.execute();

    auto res_open = sort_2.output.openRelation();

    bool test_res = checkCorrectness(passwords_data, res_open,
                                     0, 1, 3, 2,
                                     passwords_size);

    assert(test_res == true);

    // free temp data
    free(passwords_data);
    passwords.freeRelation();
}


int main(int argc, char **argv) {

    // initialize communication
    init(argc, argv);
    init_sharing();

    // Testing and timing the query
    timeTest(planner_qpwd_test, "planner_qpwd_test");

    // tear down communication
    MPI_Finalize();
    return 0;
}