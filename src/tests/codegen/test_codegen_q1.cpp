#include <assert.h>
#include <set>
#include <stdio.h>

#include "codegen-utils.h"

/***
 * Compares Secrecy's result for Q1 with result of plaintext computation.
 * @param r_data - Plaintext data for table R
 * @param s_data - Plaintext data for table S
 * @param res_open - Opened result
 * @param r_id - Index of distinct key in R
 * @param s_id - Index of join key in S
 * @param res_id - Index of R.id in `res_open`
 * @param res_sel - Index of SEL attribute in `res_open`
 * @param rows_1 - R cardinality
 * @param rows_2 - S cardinality
 * @param rows_3 - Output cardinality
 * @return True if the two results are the same, False otherwise.
 */
bool checkQ1Correctness(Data **r_data, Data **s_data,
                        RowData res_open,
                        int r_id, int s_id, int res_id, int res_sel,
                        int rows_1, int rows_2, int rows_3) {

    if (get_rank() == 0) {
        // Collect distinct keys from R
        std::set<Data> r_b_data_set;
        for (int i = 0; i < rows_1; ++i) {
            r_b_data_set.insert(r_data[i][r_id]);
        }

        // Collect distinct keys from S
        std::set<Data> s_b_data_set;
        for (int i = 0; i < rows_2; ++i) {
            s_b_data_set.insert(s_data[i][s_id]);
        }

        // Collect common keys
        std::set<Data> intersection;
        set_intersection(r_b_data_set.begin(), r_b_data_set.end(), s_b_data_set.begin(),
                         s_b_data_set.end(), std::inserter(intersection, intersection.begin()));

        // Collect distinct keys from secure execution
        auto res_open_ptr = res_open.get();
        std::set<Data> res_set;
        for (int i = 0; i < rows_3; ++i) {
            if (res_open_ptr[res_sel][i]) {
                if (res_open_ptr[res_id][i] >= 0) { // If the id is not masked
                    res_set.insert(res_open_ptr[res_id][i]);
                }
            }
        }
        printf("Distinct count %ld==%ld\n", intersection.size(), res_set.size());

        return intersection == res_set;
    } else {
        return true;
    }
}

// Generates the optimal plan and executes the query
void executeAndTestQuery(int argc, char *argv[]) {
    assert(argc==3); // Executable, left input cardinality, right input cardinality
    int rows[2] = {std::stoi(argv[1]), std::stoi(argv[2])};  // Input cardinalities
    std::string query_file = "../examples/queries/q1.txt", // Query file
                schema_file = "../examples/schemas/q1_schema.txt"; // Schema file
    std::shared_ptr<Database> db;
    std::shared_ptr<Query> query;
    // Create (empty) database and parse query
    createDBandQuery(rows, 2, query_file, schema_file, "Q1", query, "Q1_DB", db);
    // Print DB schema to stdout
    db->printSchema();

    // Populate relations with random data
    std::vector<DataTable> r_data;
    populateRandDB(db, r_data);

    // Parse query and generate plan
    query->generatePlan();
    assert(!query->plans.empty());

    // Exchange seeds
    exchange_rsz_seeds();

    // Call execution on the root operator
    leaderLogging("[INFO] Executing the plan\n");
    query->root->execute();
    leaderLogging("[INFO] Done.\n");

    // Test correctness
    auto res_open = query->root->output.openRelation();
    int out_size = query->root->output.getSize();

    auto res_open_ptr = res_open.get();

    leaderLogging("[CodeGen] Open output relation: " + std::to_string(out_size) + "\n");

    std::vector<std::shared_ptr<Relation>> schema = db->getRelations();
    const int r_id = (*schema[0])["R.[id]"].getIndex();
    const int s_id = (*schema[1])["S.[id]"].getIndex();
    const int res_id = query->root->output["R.[id]"].getIndex();
    const int res_sel = query->root->output["[SEL]"].getIndex();

    if (get_rank() == 0) {
        printf("r_id: %d, s_id: %d, res_id: %d, res_sel: %d\n", r_id, s_id, res_id, res_sel);
    }

    bool test_res = checkQ1Correctness(r_data[0], r_data[1],
                                       res_open,
                                       r_id, s_id, res_id, res_sel,
                                       rows[0], rows[1], out_size);
    assert(test_res);

    // Free input data
    free(r_data[0]);
    free(r_data[1]);
    schema[0]->freeRelation();
    schema[1]->freeRelation();
}

// Query SQL:
// SELECT DISTINCT R.id
// FROM R, S
// WHERE R.id = S.id
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