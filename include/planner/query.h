#ifndef SECRECY_QUERY_H
#define SECRECY_QUERY_H

#include "operators/operators.h"
#include "parser/parser.h"
#include "queryGraph.h"
#include "plan.h"
#include "database.h"
#include "planGen.h"
#include "codeGen.h"
#include "preprocessor.h"

/***
 * A Secrecy query.
 */
class Query {
private:
    // The database the query is applied to
    Database database;
    // The query parser
    Parser parser;
    // The query preprocessor
    Preprocessor prep;
    // The plan generator
    PlanGen pGen;
    // The code generator
    CodeGen cGen;
    // The query graph
    std::shared_ptr<QueryGraph> qGraph;

public:
    // The (optimal) logical plan(s)
    std::vector<std::shared_ptr<Plan>> plans;
    // The physical plan
    std::shared_ptr<Operator> root;

    std::vector<std::shared_ptr<Operator>> ops;
    /***
     * Query constructor.
     * @param filename
     * @param db
     * @param fileInput
     */
    Query(std::string filename, Database &db, bool fileInput = false) : parser(Parser(filename, fileInput)),
                                                                        database(db) {
        // Parse query string and do syntax checking
        this->parser.parseQuery();
    }

    /***
     * Generates optimal query plan.
     */
    void generatePlan(bool physical=true) {
        // Preprocess query and create the query graph
        leaderLogging("[INFO] Preprocessing...\n");
        this->qGraph = this->prep.transform(this->parser.query, this->database);
        leaderLogging("Done.\n");

        // Generate optimal plan
        leaderLogging("[INFO] Planning...\n");
        this->pGen.generatePlan(this->qGraph, this->plans, this->database);
        leaderLogging("[INFO] Done.\n");

        if (physical) {
            leaderLogging("[INFO] Generating physical plan...\n");
            root = this->cGen.generatePhysicalPlan(plans, &database);
            leaderLogging("[INFO] Done.\n");
        }
    }

    void printPlan() {
        for (auto plan : plans) {
            if (!plan) {
                std::cout << "No plan found." << std::endl;
                return;
            }
            plan->print();
        }
    }

};

#endif //SECRECY_QUERY_H
