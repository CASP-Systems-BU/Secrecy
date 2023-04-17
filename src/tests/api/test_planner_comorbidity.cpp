#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"
#include <map>


bool checkComorbidityCorrectness(Data **diagnosis_data, int diagnosis_size, 
                                Data **cdiff_data, int cdiff_size,
                                Relation relation) {

    // SELECT diag, COUNT(*) cnt
    // FROM diagnosis
    // WHERE pid IN cdiff_cohort
    // GROUP BY diag
    // ORDER BY cnt DESC
    // LIMIT 10
    
    auto res_open = relation.openRelation();
    auto res_open_ptr = res_open.get();

    const int d_pid = 0;
    const int d_diag = 1;

    const int c_cdiff = 0;

    const int res_diag = relation["[diag]"].getIndex();
    const int res_count = relation["[count]"].getIndex();

    if (get_rank() == 0) {
        std::map<Data,Data> count_real;
        for(int i = 0; i < diagnosis_size; ++i){
            for(int j = 0; j < cdiff_size; ++j){
                if(diagnosis_data[i][d_pid] == cdiff_data[j][c_cdiff]){
                    count_real[diagnosis_data[i][d_diag]]++;
                    break;
                }
            }
        }

        // convert to vector and sort it
        std::vector<std::pair<Data, Data> > vec_count_real;
        for(auto it = count_real.begin(); it!=count_real.end(); ++it){
            vec_count_real.push_back(std::pair<Data, Data> (it->second, it->first));
        }
        std::sort(vec_count_real.begin(), vec_count_real.end());


        std::vector<std::pair<Data, Data> > vec_count_res;
        for(int i = 0; i < relation.getSize(); ++i){
            if(res_open_ptr[res_diag][i] != -1){
                vec_count_res.push_back(std::pair<Data, Data>(res_open_ptr[res_count][i], res_open_ptr[res_diag][i]));
            }
        }

        // for (int i =0; i < vec_count_real.size(); ++i){
        //     printf("%ld==%ld\t\t%ld==%ld\n", vec_count_real[i].first, vec_count_res[i].first, vec_count_real[i].second, vec_count_res[i].second);
        // }

        printf("Groups count %ld==%ld\n", vec_count_real.size(), vec_count_res.size());

        return vec_count_real == vec_count_res;
    } else {
        return true;
    }
}

void planner_comorbidity_test (){
    // Diagnosis Table
    int diagnosis_size = 128;
    std::vector<std::string> diagnosis_schema = {
        "[pid]", "[diag]", "[count]", "[SEL]", "count", "SEL"
    };
    Data **diagnosis_data = getRandomData(diagnosis_size, diagnosis_schema.size());
    Relation diagnosis = getTestRelation("DIAGNOSIS", diagnosis_size, diagnosis_schema, diagnosis_data);

    // Patients Table
    int cdiff_size  = 128;
    std::vector<std::string> cdiff_schema = {
        "[cdiff]"
    };
    Data ** cdiff_data = getRandomData(cdiff_size, cdiff_schema.size());
    Relation cdiff = getTestRelation("CDIFF", cdiff_size, cdiff_schema, cdiff_data);

    auto uid_expr = Expression("DIAGNOSIS", "[pid]");
    auto diag_expr = Expression("DIAGNOSIS", "[diag]");
    auto count_expr_a = Expression("DIAGNOSIS", "count");
    auto count_expr_b = Expression("DIAGNOSIS", "[count]");
    auto sel_a = Expression("DIAGNOSIS", "SEL");
    auto sel_b= Expression("DIAGNOSIS", "[SEL]");

    auto cdiff_expr = Expression("CDIFF", "[cdiff]");


    //exchange seeds
    exchange_rsz_seeds();

    // SELECT diag, COUNT(*) cnt
    // FROM diagnosis
    // WHERE pid IN cdiff_cohort
    // GROUP BY diag
    // ORDER BY cnt DESC
    // LIMIT 10

    // Optimized Query
    ScanOperator diagnosis_Scan(diagnosis, 4);
    ScanOperator cdiff_Scan(cdiff, 4);
    // SortOperator sort(diagnosis_Scan, {uid_expr}, {0}, 4);
    SemiJoinOperator semiJoin(diagnosis_Scan, cdiff_Scan, uid_expr, cdiff_expr, 4);
    Expression star(ExpressionType::Star);
    GroupByOperator groupBy(semiJoin, GroupByType::COUNT, {sel_b, diag_expr}, {count_expr_a, count_expr_b}, {GroupByType::COUNT}, {star}, true, 4);
    ConversionOperator conversion(groupBy, {count_expr_a}, {count_expr_b});
    // TODO: Sorting on diag just for correctness test.
    SortOperator sort(conversion, {count_expr_b, diag_expr}, {1, 1}, 4);

    // Execution
    sort.execute();
    auto output_relation = sort.output;

    // Checking result
    bool test_res = checkComorbidityCorrectness(diagnosis_data, diagnosis_size,
                                                cdiff_data, cdiff_size,
                                                output_relation);

    assert(test_res == true);

    // Free temp data
    free(diagnosis_data);
    free(cdiff_data);
    diagnosis.freeRelation();
    cdiff.freeRelation();
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