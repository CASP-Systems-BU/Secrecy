#include <assert.h>
#include <stdio.h>

#include "../test-utils.h"


// Parses program arguments and creates database and query objects
static void createDBandQuery(int argc, char *argv[], const std::string &db_name, std::shared_ptr<Database> &db,
                             std::shared_ptr<Query> &q) {
    assert(argc>=4);  // 1 executable, 1 query filename, 1 schema filename, and at least 1 input relation size
    // Parse arguments
    std::string query_file, schema_file;
    int rel_no = argc-3;
    int rows[rel_no];
    parseArguments(argc, argv, query_file, schema_file, rows, rel_no);
    // Parse database schema
    std::vector<std::shared_ptr<Relation>> schema;
    create_schema(schema_file, rows, rel_no, schema);
    // Create (empty) database
    db = std::make_shared<Database>(Database(db_name, schema));
    // Parse query
    q = std::make_shared<Query>(Query(query_file, *db, true));
}

// Creates database and query objects for the given query file and database schema
static void createDBandQuery(int *rows, int rel_no,
                             const std::string &query_file, const std::string &schema_file,
                             const std::string &q_name, std::shared_ptr<Query> &q,
                             const std::string &db_name, std::shared_ptr<Database> &db) {
    // Create (empty) database and parse query
    leaderLogging("[INFO] Creating DB and parsing query "+ q_name +" ...\n");
    std::vector<std::shared_ptr<Relation>> schema;
    create_schema(schema_file, rows, rel_no, schema);
    db = std::make_shared<Database>(Database(db_name, schema));
    q = std::make_shared<Query>(Query(query_file, *db, true));
    leaderLogging("[INFO] Done.\n");
}

// Populates tables in the given database with random secret-shared data
static void populateRandDB(std::shared_ptr<Database> db, std::vector<DataTable> &r_data) {
    std::vector<std::shared_ptr<Relation>> relations = db->getRelations();
    int rel_no = relations.size();
    r_data.resize(rel_no);
    for (int i=0; i<rel_no; i++) {
        // Generate random data
        r_data[i] = getRandomData(relations[i]->getSize(), relations[i]->colsNum());
        // Secret share random data and populate relations
        populateTestRelation(relations[i], r_data[i]);
    }
}
