#include "../../include/benchmarking/data.h"

#include "../../include/common/communication.h"
#include <sys/time.h>
#include <iostream>
#include <fstream>

Data **getRandomData(int rows, int cols) {
    Data **data = allocateContent<Data>(rows, cols, 1);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            data[i][j] = random() % 100;
        }
    }
    return data;
}



Data ** readDataFromFile(const std::string& filename, const int& rows, const int& cols) {

    Data **data = allocateContent<Data>(rows, cols, 1);
    std::ifstream file(filename);
    if (file.is_open()) {
        int num_cols = 0;
        char c;
        
        file >> c;
        if (c != '{') throw std::runtime_error("Invalid format");
        for (int i = 0; i < rows; i++) { // loop over rows
            file >> c;
            if (c != '{') throw std::runtime_error("Invalid format");
            for (int j = 0; j < cols; j++) { // loop over columns
                file >> data[i][j];
                file >> c;
                if (j == cols-1 && c != '}') throw std::runtime_error("Invalid format");
                else if (j < cols-1 && c != ',') throw std::runtime_error("Invalid format");
            }
            file >> c;
            if (i < rows-1 && c != ',') throw std::runtime_error("Invalid format");
        }
        if (c != '}') throw std::runtime_error("Invalid format");
    } else {
        throw std::runtime_error("Unable to open file");
    }
    return data;
}


// Same as getTestRelation but populates an existing Relation
// with test data
void populateTestRelation(std::shared_ptr<Relation> &relation, Data **data) {
    int rows = relation->getSize();
    auto cols = relation->getCols(true);
    int num_cols = cols.size();

    // AShare part
    Share **shares_1 = allocateContent<AShare>(rows, num_cols);
    if (get_rank() == 0) {
        // Secret share the dataA
        Share **shares_2 = allocateContent<AShare>(rows, num_cols);
        Share **shares_3 = allocateContent<AShare>(rows, num_cols);

        Share element, s1, s2, s3;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < num_cols; ++j) {
                element = data[i][j];

                if (relation->bShares[j]) {
                    generate_bool_share(element, &s1, &s2, &s3);
                } else {
                    generate_int_share(element, &s1, &s2, &s3);
                }

                int col_1 = R_3PC.PARTY_REPLICATION * j;
                int col_2 = col_1 + 1;

                // 1st table:
                shares_1[i][col_1] = s1;
                shares_1[i][col_2] = s2;

                // 2nd table:
                shares_2[i][col_1] = s2;
                shares_2[i][col_2] = s3;

                // 3rd table:
                shares_3[i][col_1] = s3;
                shares_3[i][col_2] = s1;
            }
        }
        //  Send table to other parties
        sendArrayToParty<AShare>(shares_2[0], rows, cols.size(), 1, TEST_SHARE_TAG);
        sendArrayToParty<AShare>(shares_3[0], rows, cols.size(), 2, TEST_SHARE_TAG);
        free(shares_2);
        free(shares_3);
    } else {
        // Receive table from the leader
        getArrayFromParty<AShare>(shares_1[0], rows, cols.size(), 0, TEST_SHARE_TAG);
    }
    AShareTable shareTable = newShareTable<Share, ShareTable>(rows, cols.size(), shares_1);
    relation->insertTable(shareTable);
    freeShareTable<AShareTable>(shareTable);
}

Relation getTestRelation(std::string name, int rows,
                         std::vector<std::string> cols, Data **data) {

    Relation relation = Relation(name, cols, rows);

    // AShare part
    Share **shares_1 = allocateContent<AShare>(rows, cols.size());
    if (get_rank() == 0) {
        // Secret share the dataA
        Share **shares_2 = allocateContent<AShare>(rows, cols.size());
        Share **shares_3 = allocateContent<AShare>(rows, cols.size());

        Share element, s1, s2, s3;
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols.size(); ++j) {
                element = data[i][j];

                // TODO: create table of functions pointer
                //  to assign share generation function for each column
                //  just once at the beginning.
                if (relation.isBColumnName(cols[j])) {
                    generate_bool_share(element, &s1, &s2, &s3);
                } else {
                    generate_int_share(element, &s1, &s2, &s3);
                }

                int col_1 = R_3PC.PARTY_REPLICATION * j;
                int col_2 = col_1 + 1;

                // 1st table:
                shares_1[i][col_1] = s1;
                shares_1[i][col_2] = s2;

                // 2nd table:
                shares_2[i][col_1] = s2;
                shares_2[i][col_2] = s3;

                // 3rd table:
                shares_3[i][col_1] = s3;
                shares_3[i][col_2] = s1;
            }
        }
        //  Send table to other parties
        sendArrayToParty<AShare>(shares_2[0], rows, cols.size(), 1, TEST_SHARE_TAG);
        sendArrayToParty<AShare>(shares_3[0], rows, cols.size(), 2, TEST_SHARE_TAG);
        free(shares_2);
        free(shares_3);
    } else {
        // Receive table from the leader
        getArrayFromParty<AShare>(shares_1[0], rows, cols.size(), 0, TEST_SHARE_TAG);
    }
    AShareTable shareTable = newShareTable<Share, ShareTable>(rows, cols.size(), shares_1);
    relation.insertTable(shareTable);
    freeShareTable<AShareTable>(shareTable);

    return relation;
}


void timeTest(void (func)(), std::string name) {
    struct timeval begin, end;
    long seconds, micro;
    double elapsed;

    // start timer
    gettimeofday(&begin, 0);
    func();

    // stop timer
    gettimeofday(&end, 0);
    seconds = end.tv_sec - begin.tv_sec;
    micro = end.tv_usec - begin.tv_usec;
    elapsed = seconds + micro * 1e-6;

    printf("%s \tpassed in:\t\t%f\n", name.c_str(), elapsed);
}
