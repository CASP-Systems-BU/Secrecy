#include <assert.h>
#include <map>
#include <stdio.h>

#include "codegen-utils.h"

/***
 * Compares Secrecy's result for Q2 with result of plaintext computation.
 * @param r_data
 * @param s_data
 * @param res_open
 * @param r_id
 * @param r_ak - Index of attribute ak in R
 * @param s_id
 * @param res_ak - Index of attribute ak in `res_open`
 * @param res_sel
 * @param res_sel_a - Index of the arithmetic SEL attribute in `res_open`
 * @param rows_1
 * @param rows_2
 * @param rows_3
 * @return
 */
bool checkQ2Correctness(Data **r_data, Data **s_data,
                        RowData res_open,
                        int r_id, int r_ak, int s_id,
                        int res_ak, int res_sel, int res_sel_a,
                        int rows_1, int rows_2, int rows_3) {

    if (get_rank() == 0) {

        // Count the number of join matches per distinct r_ak attribute (only for those ak value that have a match)
        std::map<Data, long long> real_count;
        for (int i_1 = 0; i_1 < rows_1; ++i_1) {
            for (int i_2 = 0; i_2 < rows_2; ++i_2) {
                if (r_data[i_1][r_id] == s_data[i_2][s_id]) {
                    real_count[r_data[i_1][r_ak]]++;
                }
            }
        }

        // Collect the result from the secure execution
        auto res_open_ptr = res_open.get();
        std::map<Data, long long> res_count;
        for (int i = 0; i < rows_3; ++i) {
            if (res_open_ptr[res_ak][i] >= 0) {
                res_count[res_open_ptr[res_ak][i]] = res_open_ptr[res_sel_a][i];
            }
        }

        // Add r_ak with no match with a zero value
        for(auto it = res_count.cbegin(); it != res_count.cend(); ++it)
            real_count[it->first];

        printf("Groups count %ld==%ld\n", res_count.size(), real_count.size());
        return real_count == res_count;
    } else {
        return true;
    }
}

// Generates the optimal plan and executes the query
void executeAndTestQuery(int argc, char *argv[]) {
    assert(argc==3); // Executable, left input cardinality, right input cardinality
    int rows[2] = {std::stoi(argv[1]), std::stoi(argv[2])};  // Input cardinalities
    std::string query_file = "../examples/queries/q2.txt", // Query file
    schema_file = "../examples/schemas/q2_schema.txt"; // Schema file
    std::shared_ptr<Database> db;
    std::shared_ptr<Query> query;
    // Create (empty) database and parse query
    createDBandQuery(rows, 2, query_file, schema_file, "Q2", query, "Q2_DB", db);
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
    leaderLogging("[CodeGen] Executing the plan\n");
    query->root->execute();
    leaderLogging("[INFO] Done.\n");

    // Open and test correctness
    auto res_open = query->root->output.openRelation();
    int out_size = query->root->output.getSize();
    auto res_open_ptr = res_open.get();

    std::vector<std::shared_ptr<Relation>> schema = db->getRelations();
    const int r_id = (*schema[0])["R.[id]"].getIndex();
    const int r_ak = (*schema[0])["R.[ak]"].getIndex();
    const int s_id = (*schema[1])["S.[id]"].getIndex();
    const int res_ak = query->root->output["R.[ak]"].getIndex();
    const int res_sel = query->root->output["R.[SEL]"].getIndex();
    const int res_sel_a = query->root->output["R.SEL"].getIndex();

    bool test_res = checkQ2Correctness(r_data[0], r_data[1],
                                       res_open,
                                       r_id, r_ak, s_id,
                                       res_ak, res_sel, res_sel_a,
                                       rows[0], rows[1], out_size);
    assert(test_res);

    // Free input data
    free(r_data[0]);
    free(r_data[1]);
    schema[0]->freeRelation();
    schema[1]->freeRelation();
}

// SELECT R.ak , COUNT(*)
// FROM R, S
// WHERE R.id = S.id
// GROUP BY R.ak
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