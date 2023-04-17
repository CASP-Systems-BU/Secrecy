#ifndef SECRECY_DATA_H
#define SECRECY_DATA_H

#include "../planner/relation.h"
#include <vector>

// TODO: add function that gets random data by number of rows and
//  list of columns names so it counts for "COL" and "[COL]"

Data **getRandomData(int rows, int cols);

/**
 * Read 2D data from file and returns a 2D array accessible as [row_index][col_index].
 * Acceptable file format is like "{{1,2,3}, {4,5,6}, {7,8,9}}"; spaces and new lines
 * do not affect correctness of parsing.
 *  Usage example auto data = readDataFromFile("../data/q2/r_1.txt", 3 , 3); data[i][j] = 0;
 * @param filename path to file with the 2D data.
 * @param rows number of rows in the 2D data.
 * @param cols number of cols in each row in the 2D data. Must be equal between all rows or an error is thrown.
 * @return 2D array accessible as [row_index][col_index]
 */
Data **readDataFromFile(const std::string& filename, const int& rows, const int& cols);


// get test relations
Relation getTestRelation(std::string name, int rows,
                         std::vector<std::string> cols, Data **data);

void populateTestRelation(std::shared_ptr<Relation> &relation, Data **data);

void timeTest(void (func)(), std::string name);

#endif //SECRECY_DATA_H