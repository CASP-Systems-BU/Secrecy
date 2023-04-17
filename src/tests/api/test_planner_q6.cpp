#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"
#include <set>

// SELECT SUM(L_EXTENDEDPRICE*L_DISCOUNT) AS REVENUE
// FROM LINEITEM
// WHERE L_SHIPDATE >= '1994-01-01'         == 1
// AND L_SHIPDATE < dateadd(yy, 1, cast('1994-01-01' as date))          == 100
// AND L_DISCOUNT BETWEEN .06 - 0.01 AND .06 + 0.01
// AND L_QUANTITY < 24

// TODO: Update the plain text evaluation

/**
 *
 * @param line_data Plain text data for the table "line"
 * @param rows_2 Size for table "line"
 * @param line_schema Schema for table "line"
 * @param res Data return as output
 * @return
 */
bool checkQ6Correctness(Data **line_data, const int & rows_2,
                        const std::vector<std::string>& line_schema,
                        Data res) {

    if (get_rank() == 0) {
        // 0, 1, 7, 8 
        const int l_ship_date = 0;
        const int l_quantity = 1;
        const int l_discount = 7;
        const int l_extend_price = 8;

        Data res_real = 0;

        for(int i =0 ; i < rows_2; ++i){
            if (line_data[i][l_ship_date] >= 1 && line_data[i][l_ship_date] < 100
                && line_data[i][l_discount] >= 1 && line_data[i][l_discount] <= 6
                && line_data[i][l_quantity] < 24){
                    res_real+= (line_data[i][l_discount]* line_data[i][l_extend_price]);
                }
        }


        // First table set
        printf("%lld==%lld\n", res_real, res);
        return res == res_real;
    } else {
        return true;
    }
}


// NOTE: two ways to implement this query
// - calculate ashare("[SEL]") .. ("SEL"*"extend_price"*"discount").sum()
// - group by sum on "[SEL]" and ("extend_price"*"discount")
// TODO: need to calculate .sum()
// TODO: need to track the selection bit:
//  - track size: "long long" or "byte"
//  - create function to keep track that [SEL] and SEL are synced
//  - update the to be using the set size for selection bit.
//  - Better to have a conversion function that has number of bits as input parameters.
void planner_q6_test() {

    // Line table
    int line_row = 128;
    std::vector<std::string> line_schema = {"[ship_date]", "[quantity]",
                                            "[0-ship_date]", "[ship_date-100]",
                                            "[quantity-24]",
                                            "[discount-1]", "[6-discount]",
                                            "discount", "extend_price"};
    Data **line_data = getRandomData(line_row, line_schema.size());
    for(int i =0; i < line_row; ++i){
        line_data[i][2] = 0 - line_data[i][0];
        line_data[i][3] = line_data[i][0] - 100;
        line_data[i][4] = line_data[i][1] - 24;
        line_data[i][5] = line_data[i][7] - 1;
        line_data[i][6] = 6 - line_data[i][7];
    }
    Relation line = getTestRelation("LINE", line_row, line_schema, line_data);

    //exchange seeds
    exchange_rsz_seeds();

    // QUERY OPTIMIZATION:
    Column& sel_b = line["[0-ship_date]"] < 0 & line["[ship_date-100]"] < 0
            & line["[quantity-24]"] < 0 & line["[discount-1]"] >= 0
            & line["[6-discount]"] >= 0;
    AColumn<AShare, AShareTable> sel_a  = convertToAColumn<AShare, AShareTable>(sel_b)
            * line["discount"] * line["extend_price"];
    
    AShare result = sel_a.sum();


    Data opened_result = open_a(result);

    bool test_res = checkQ6Correctness(line_data, line_row,
                                        line_schema, opened_result);

    assert(test_res == true);

    // free temp data
    free(line_data);
    line.freeRelation();
}


int main(int argc, char **argv) {

    // initialize communication
    init(argc, argv);
    init_sharing();

    // Testing and timing the query
    timeTest(planner_q6_test, "planner_q6_test");

    // tear down communication
    MPI_Finalize();
    return 0;
}