#ifndef SECRECY_PARSER_H
#define SECRECY_PARSER_H

#include <fstream>
#include <streambuf>
#include <iostream>

#include "../../external-lib/sql-parser/src/SQLParser.h"
#include "../../external-lib/sql-parser/src/util/sqlhelper.h"

/***
 * A wrapper to the Hyrise SQL parser.
 */
class Parser {
private:
    hsql::SQLParserResult result;
public:
    /***
     * A parsed SQL query.
     */
    struct SQLQuery {
        /**
         * The query id
         */
        int id;
        /**
         * The query string in SQL syntax.
         */
        std::string qString;
        /**
         * A pointer to the Abstract Syntax Tree (AST).
         */
        const hsql::SelectStatement *tree;
        /**
         * Pointers to ASTs of nested queries.
         */
        std::vector<const hsql::SelectStatement*> subTrees;
        /**
         * Returns a pointer to the SELECT clause of the query (if any).
         * @return A pointer to the AST of the SELECT clause
         */
        const std::vector<hsql::Expr*>* selectClause() const { return this->tree->selectList; };
        /**
         * @return True if there is a SELECT DISTINCT clause, False otherwise.
         */
        bool hasDistinct() const { return this->tree->selectDistinct; };
        /**
         * Returns a pointer to the FROM clause of the query.
         * @return A pointer to the AST of the FROM clause
         */
        const hsql::TableRef* fromClause() const { return this->tree->fromTable; };
        /**
         * Returns a pointer to the WHERE clause of the query.
         * @return A pointer to the AST of the WHERE clause
         */
        const hsql::Expr* whereClause() const { return this->tree->whereClause; };
        /**
         * Returns a pointer to the GROUP-BY clause of the query.
         * @return A pointer to the AST of the GROUP-BY clause
         */
        const hsql::GroupByDescription* groupByClause() const { return this->tree->groupBy; };
        /**
         * Returns a pointer to the ORDER-BY clause of the query.
         * @return A pointer to the AST of the ORDER-BY clause
         */
        const std::vector<hsql::OrderDescription*>* orderByClause() const { return this->tree->order; };
        /**
         * Returns a pointer to the WITH clause of the query.
         * @return A pointer to the AST of the WITH clause
         */
        const std::vector<hsql::WithDescription*>* withClause() const { return this->tree->withDescriptions; };
    };

    /***
     * The parsed query.
     */
    SQLQuery query;
    /***
     * The query file name.
     */
    std::string fileName;
    /***
     * True if parsing was successful, False if no query has been parsed yet.
     */
    bool parsingDone;

    /***
     * Constructor.
     * @param input - The query string or filename
     * @param fileInput - True if `input` is a filename, False if `input` is the query string itself.
     */
    Parser(std::string input, bool fileInput = false) {
        if (fileInput) {
            this->fileName = input;
            std::ifstream inputStream(input);

            if (!inputStream) {
                leaderLogging("[ERROR] Could not open file '" + input + "'\n", true);
            }

            this->query.qString = std::string((std::istreambuf_iterator<char>(inputStream)),
                                              std::istreambuf_iterator<char>());
        } else {
            this->query.qString = input;
        }
        // No
        this->parsingDone = false;
    }

    bool parseQuery() {
        if (!this->parsingDone) {
            bool res = hsql::SQLParser::parse(this->query.qString, &result);
            if (!this->result.isValid()) {
                leaderLogging("[ERROR] Given string is not a valid SQL query.\n");
                leaderLogging(this->result.errorMsg() +
                                this->result.errorLine() +
                                this->result.errorColumn(), true);
            }
            auto q = this->result.getStatement(0);
            this->query.tree = (hsql::SelectStatement *) q;
            this->parsingDone = true;
        }
        return true;
    }

    void printQueryRepresentation() {

        if (result.isValid()) {
            printf("Parsed successfully!\n");
            printf("Number of statements: %lu\n", result.size());

            for (auto i = 0u; i < result.size(); ++i) {
                // Print a statement summary.
                hsql::printStatementInfo(result.getStatement(i));
            }
        } else {
            fprintf(stderr, "Given string is not a valid SQL query.\n");
            fprintf(stderr, "%s (L%d:%d)\n",
                    result.errorMsg(),
                    result.errorLine(),
                    result.errorColumn());
        }
    }

    void printDebugging() {

        auto s1 = this->result.getStatement(0);
        printStatementInfo(s1);

        hsql::SelectStatement *stmt = (hsql::SelectStatement *) s1;
        if (stmt->selectDistinct) {
            std::cout << "We have distinct:" << std::endl;
            for (auto expr : stmt->selectList[0]) {
                printExpression(expr, 4);
                if (expr->distinct) {
                    std::cout << "this is distinct \n";
                }
            }
        }

        if (stmt->limit != nullptr) {
            if (stmt->limit->limit != nullptr) {
                int limit = stmt->limit->limit->ival;
                std::cout << limit << std::endl;
            } else {
                std::cout << "no limit" << std::endl;
            }
        } else {
            std::cout << "no limit" << std::endl;
        }

        if (stmt->fromTable != nullptr) {
            if (stmt->fromTable->type == hsql::kTableJoin) {
                std::cout << ("We have a join here.\n");
            }

            if (stmt->fromTable->name != nullptr) {
                std::cout << ("One table with name") << std::endl;
            }
        } else {
            std::cout << ("from table is null") << std::endl;
        }
    }
};

#endif //SECRECY_PARSER_H
