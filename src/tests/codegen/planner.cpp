#include "codegen-utils.h"

// Generates the optimal plan for the given query
void executeQuery(int argc, char *argv[]){
    std::shared_ptr<Database> db;
    std::shared_ptr<Query> query;
    // Create (empty) database and parse query
    leaderLogging("[INFO] Creating DB and parsing query...\n");
    createDBandQuery(argc, argv, "Example_DB", db, query);
    leaderLogging("[INFO] Done.\n");
    // Print DB schema to stdout
    db->printSchema();

    // Populate database relations with random data
    std::vector<DataTable> r_data;
    populateRandDB(db, r_data);

    // Parse query and generate plan
    query->generatePlan(false);
    assert(!query->plans.empty());

    // Exchange seeds
    exchange_rsz_seeds();

    // Free random input data
    auto schema = db->getRelations();
    int rel_no = schema.size();
    for (int i=0; i<rel_no; i++) {
        free(r_data[i]);
        schema[i]->freeRelation();
    }
}


int main(int argc, char *argv[]) {
    // Initialize MPI communication
    init(argc, argv);
    init_sharing();
    // Compile and execute the query
    executeQuery(argc, argv);
    // Shut down
    MPI_Finalize();
    return 0;
}