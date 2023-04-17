#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"
#include <set>

// TODO:
//  1 - Increase size
//  2 - Support column aggregations

bool checkCorrectness(Data **diagnosis_data, int diagnosis_size,
                                 Data **medications_data, int medications_size,
                                 Relation relation) {

    auto res_open = relation.openRelation();
    auto res_open_ptr = res_open.get();

    const int d_pid = 0;
    const int d_time = 1;
    const int d_diag = 2;

    const int m_pid = 0;
    const int m_time = 1;
    const int m_med = 2;

    auto pid_expr_1 = Expression("DIAGNOSIS", "[pid]");
    const int res_pid = relation[pid_expr_1].getIndex();
    const int res_sel = relation["[SEL]"].getIndex();


    if (get_rank() == 0) {

        std::set<Data> real_count;
        for(int i = 0; i < diagnosis_size; ++i){
            for(int j = 0; j < medications_size; ++j){
                if(diagnosis_data[i][d_pid] == medications_data[j][m_pid]
                && diagnosis_data[i][d_diag] == 50
                && medications_data[j][m_med] == 50
                && diagnosis_data[i][d_time] <= medications_data[j][m_time]){
                    real_count.insert(diagnosis_data[i][d_pid]);
                }
            }
        }


        std::set<Data> res_count;
        for(int i = 0; i < relation.getSize(); ++i){
            if(res_open_ptr[res_sel][i]){
                res_count.insert(res_open_ptr[res_pid][i]);
            }
        }

        printf("Distinct count %ld==%ld\n", real_count.size(), res_count.size());

        return real_count == res_count;
    } else {
        return true;
    }
}

void planner_comorbidity_test (){
    // Diagnosis Table
    int diagnosis_size = 128;
    std::vector<std::string> diagnosis_schema = {
            "[pid]", "[time]",
            "[diag]", "[diag-c]", "[c-diag]",
            "[SEL]"
    };
    Data **diagnosis_data = getRandomData(diagnosis_size, diagnosis_schema.size());
    for (int i = 0; i < diagnosis_size; ++i) {
        diagnosis_data[i][3] = diagnosis_data[i][2] - 50;
        diagnosis_data[i][4] = 50 - diagnosis_data[i][2];
    }
    Relation diagnosis = getTestRelation("DIAGNOSIS", diagnosis_size, diagnosis_schema, diagnosis_data);

    // Patients Table
    int medications_size  = 128;
    std::vector<std::string> medications_schema = {
            "[pid]", "[time]",
            "[med]", "[med-c]", "[c-med]",
            "[SEL]", "SEL"
    };
    Data ** medications_data = getRandomData(medications_size, medications_schema.size());
    for (int i = 0; i < diagnosis_size; ++i) {
        medications_data[i][3] = medications_data[i][2] - 50;
        medications_data[i][4] = 50 - medications_data[i][2];
    }
    Relation medications = getTestRelation("MEDICATIONS", medications_size, medications_schema, medications_data);

    auto pid_expr_1 = Expression("DIAGNOSIS", "[pid]");
    auto pid_expr_2 = Expression("MEDICATIONS", "[pid]");
    auto time_expr_1 = Expression("DIAGNOSIS", "[time]");
    auto time_expr_2 = Expression("MEDICATIONS", "[time]");
    auto sel_b_1 = Expression("DIAGNOSIS", "[SEL]");
    auto sel_b_2 = Expression("MEDICATIONS", "[SEL]");
    auto sel_b_3 = Expression("DIAGNOSIS*MEDICATIONS", "[SEL]");
    auto sel_a_3 = Expression("DIAGNOSIS*MEDICATIONS", "SEL");

    auto join_expr = Expression(ExpressionType::Equal, pid_expr_1, pid_expr_2);

    //exchange seeds
    exchange_rsz_seeds();

    // SELECT count(DISTINCT pid)
    // FROM diagnosis as d, medication as m on d.pid = m.pid
    // WHERE d.diag = hd AND m.med = aspirin
    // AND d.time <= m.time

    diagnosis["[SEL]"] = !((diagnosis["[diag-c]"] < 0) ^ (diagnosis["[c-diag]"] < 0));
    diagnosis.bSelectionUsed = true;
    medications["[SEL]"] = !((medications["[med-c]"] < 0) ^ (medications["[c-med]"] < 0));
    medications.bSelectionUsed = true;

    // Optimized Query
    ScanOperator diagnosis_Scan(diagnosis, 4);
    ScanOperator medications_Scan(medications, 4);
    SortOperator sort_1(diagnosis_Scan, {sel_b_1, pid_expr_1, time_expr_1}, {1, 1, 1}, 4);
    SortOperator sort_2(medications_Scan, {sel_b_2, pid_expr_2, time_expr_2}, {1, 1, 0}, 4);
    DistinctOperator distinct_1(sort_1, {pid_expr_1}, false, 4);
    DistinctOperator distinct_2(sort_2, {pid_expr_2}, false, 4);
    JoinOperator join(distinct_1, distinct_2, join_expr, 1, 8);

    join.execute();
    auto temp_relation = join.output;
    temp_relation[sel_b_3] = temp_relation[sel_b_3] & (temp_relation[time_expr_1] <= temp_relation[time_expr_2]);

    ScanOperator scan_3(temp_relation, 4);
    ConversionOperator conversion(scan_3, {sel_b_3}, {sel_a_3});
    conversion.execute();

    auto output_relation = conversion.output;

    // Checking result
    bool test_res = checkCorrectness(diagnosis_data, diagnosis_size,
                                        medications_data, medications_size,
                                        output_relation);

    assert(test_res == true);

    // Free temp data
    free(diagnosis_data);
    free(medications_data);
    diagnosis.freeRelation();
    medications.freeRelation();
}


int main(int argc, char **argv) {

    // initialize communication
    init(argc, argv);
    init_sharing();

    // Testing and timing the query
    timeTest(planner_comorbidity_test, "planner_comorbidity_test");

    // tear down communication
    MPI_Finalize();
    return 0;
}