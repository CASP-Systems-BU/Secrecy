#include <sstream>

#include "../../include/benchmarking.h"
#include "../../include/secrecy.h"


// Trims left
static std::string ltrim(const std::string &s, const std::string& chars) {
    size_t start = s.find_first_not_of(chars);
    return (start == std::string::npos) ? "" : s.substr(start);
}

// Trims right
static std::string rtrim(const std::string &s, const std::string& chars) {
    size_t end = s.find_last_not_of(chars);
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

// Trims left and right
static std::string trim(const std::string &s, const std::string& chars) {
    return rtrim(ltrim(s, chars), chars);
}

// Extracts column name from table declaration
static std::string trimColName(const std::string& col) {
    const std::string space = " \n\r\t\f\v";
    std::string c = trim(col, space);
    c = ltrim(c, "(");
    c = rtrim(c, ")");
    c = trim(c, space);
    return c;
}

// Tokenizes a given string using `dlm` as delimiter
static void tokenize(const std::string& s, const char& dlm, std::vector<std::string>& tokens) {
    std::istringstream ss(s);
    std::string col;
    while ( getline(ss, col, dlm))
        tokens.push_back(trimColName(col));
}

// Parses program arguments
static void parseArguments(int argc, char *argv[], std::string &query_file, std::string &schema_file,
                           int *rows, int rel_no) {
    // Parse query filepath
    query_file = argv[1];
    // Parse schema filepath
    schema_file = argv[2];
    // Parse input relation sizes (in number of rows)
    for (int i = 3, j=0; i < argc; i++, j++)
        rows[j] = std::stoi(argv[i]);
    // Print arguments to stdout
    leaderLogging("[INFO] Query file: " + query_file + "'\n");
    leaderLogging("[INFO] Schema file: " + schema_file + "'\n");
    leaderLogging("[INFO] Input rows: ", rows, rel_no);
}

// Parses database schema from input file and creates (empty) relations
static void create_schema(std::string filename, int* rows, int rel_no,
                          std::vector<std::shared_ptr<Relation>>& relations) {
    std::ifstream infile(filename);
    if (!infile) {
        leaderLogging("[ERROR] Could not open file '" + filename + "'\n", true);
    }
    const std::string space = " \n\r\t\f\v";
    std::string line, rel, cols;
    std::vector<std::string> tokens, sTypes;
    std::vector<bool> b_sharing;
    int i = 0;
    while( getline(infile, line) ) {
        auto trimmed = trim(line, space);
        if (trimmed[0]=='#') { // Metadata
            sTypes.clear();
            b_sharing.clear();
            line = trimmed.substr(1);
            tokenize(line, ',', sTypes);
            b_sharing.resize(sTypes.size(), false);
            for (int k=0; k<sTypes.size(); k++)
                if (sTypes[k] == "B") b_sharing[k] = true;
        }
        else {
            tokenize(line, '(', tokens);
            rel = tokens[0];  // Relation name
            cols = tokens[1]; // Column names
            tokens.clear();
            // Extract column names
            tokenize(cols, ',', tokens);
            assert(i < rel_no);
            // Add extra SEL attributes
            tokens.push_back("SEL");
            b_sharing.push_back(true);
            // Create relation
            auto r = std::make_shared<Relation>(rel, tokens, rows[i++],
                                                                b_sharing);
            relations.push_back(r);
            tokens.clear();
        }
    }
    assert(i == rel_no); // Make sure the number of relations equals the number of cardinalities
}