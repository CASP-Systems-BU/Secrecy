#include <assert.h>
#include <set>
#include <stdio.h>

#include "codegen-utils.h"

/***
 * Compares Secrecy's result for Q3 with result of plaintext computation.
 * @param r_data
 * @param s_data
 * @param res_open
 * @param r_ak
 * @param s_ab - The index of the column in R that contains the constant
 * @param r_id
 * @param res_id
 * @param res_sel
 * @param rows_1
 * @param rows_2
 * @param rows_3
 * @return
 */
bool checkQ3Correctness(Data **r_data, const int& ab,
                        RowData res_open,
                        int r_ak, int r_id,
                        int res_id, int res_sel,
                        int rows_1, int rows_2) {

    if (get_rank() == 0) {

        // Collect distinct ids in R that pass the filter
        std::set<Data> real_set;
        for (int i_1 = 0; i_1 < rows_1; ++i_1) {
                if (r_data[i_1][r_ak] == ab) {
                    real_set.insert(r_data[i_1][r_id]);
                }
        }

        // Collect distinct ids from secure execution
        auto res_open_ptr = res_open.get();
        std::set<Data> res_set;
        for (int i = 0; i < rows_2; ++i) {
            if (res_open_ptr[res_sel][i] != 0) { // If selected
                res_set.insert(res_open_ptr[res_id][i]);
            }
        }

        printf("Real:\n");
        for(auto& ele : real_set){
            std::cout << ele << "\t";
        }
        std::cout << std::endl;

        printf("Calculated:\n");
        for(auto& ele : res_set){
            std::cout << ele << "\t";
        }
        std::cout << std::endl;

        printf("Distinct count %ld==%ld\n", real_set.size(), res_set.size());

        return real_set == res_set;
    } else {
        return true;
    }
}

// Generates the optimal plan and executes the query
void executeAndTestQuery(int argc, char *argv[]) {
    assert(argc==2); // Executable, input cardinality
    int rows[1] = {std::stoi(argv[1])};  // Input cardinalities
    std::string query_file = "../examples/queries/q3.txt", // Query file
    schema_file = "../examples/schemas/q3_schema.txt"; // Schema file
    std::shared_ptr<Database> db;
    std::shared_ptr<Query> query;
    // Create (empty) database and parse query
    createDBandQuery(rows, 1, query_file, schema_file, "Q3", query, "Q3_DB", db);
    // Print DB schema to stdout
    db->printSchema();

    // Populate relations with random data
    std::vector<DataTable> r_data;
    populateRandDB(db, r_data);

    db->printSchema();
    // Parse query and generate plan
    query->generatePlan();
    assert(!query->plans.empty());

    // Exchange seeds
    exchange_rsz_seeds();

    // Call execution on the root operator
    leaderLogging("[CodeGen] Executing the plan\n");
    query->root->execute();
    leaderLogging("[INFO] Done.\n");

    // Test correctness
    auto res_open = query->root->output.openRelation();
    int out_size = query->root->output.getSize();
    auto res_open_ptr = res_open.get();

    leaderLogging("[INFO] Open output relation: " + std::to_string(out_size) + "\n");

    std::vector<std::shared_ptr<Relation>> schema = db->getRelations();
    const int r_id = (*schema[0])["R.[id]"].getIndex();
    const int r_ak = (*schema[0])["R.[ak]"].getIndex();
    const int res_id = query->root->output["R.[id]"].getIndex();
    const int res_ak = query->root->output["R.[ak]"].getIndex();
    const int res_sel = query->root->output["R.[SEL]"].getIndex();

    const int ab = 5;

    auto correct  = checkQ3Correctness(r_data[0], ab,
                                       res_open,
                                       r_ak, r_id,
                                       res_id, res_sel,
                                       rows[0], rows[0]);

    assert(correct);

    free(r_data[0]);
    schema[0]->freeRelation();
}

// Query:
// SELECT DISTINCT R.id
// FROM R
// WHERE R.ak = 5
int main(int argc, char **argv) {
    // Initialize MPI communication
    init(argc, argv);
    init_sharing();
    // Compile and execute the query
    executeAndTestQuery(argc, argv);
    // Shut down
    MPI_Finalize();
    return 0;
}