#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"
#include <set>

/**
 * Tests the correctness of TPC-H Q4.
 *
 * SELECT O_ORDERPRIORITY, COUNT(*) AS ORDER_COUNT FROM ORDERS
 * WHERE O_ORDERDATE >= '1993-07-01'
 * AND O_ORDERDATE < dateadd(mm,3, cast('1993-07-01' as date))
 * AND EXISTS (
 *    SELECT * FROM LINEITEM WHERE L_ORDERKEY = O_ORDERKEY AND L_COMMITDATE < L_RECEIPTDATE)
 * GROUP BY O_ORDERPRIORITY
 * ORDER BY O_ORDERPRIORITY
 **/


// TODO: test Order
// TODO: updated keys indexes
/***
 *
 * @param orders_data Plain text data for table "orders"
 * @param rows_1 Size for table "orders"
 * @param orders_schema Schema for table "orders"
 * @param line_data Plain text data for table "line"
 * @param rows_2 Size for table "line"
 * @param line_schema Schema for table "line"
 * @param output Secret-Shared Output relation after execution.
 * @return
 */
bool checkQ4Correctness(Data **orders_data, const int & rows_1,
                        const std::vector<std::string>& orders_schema,
                        Data **line_data, const int & rows_2,
                        const std::vector<std::string>& line_schema,
                        const Relation & output) {
    
    auto res_open = output.openRelation();
    auto res_open_ptr = res_open.get();

    if (get_rank() == 0) {
        const int o_key = 0;
        const int o_date = 1;
        const int o_priority = 2;
        const int o_sel = 6;

        const int l_key = 0;
        const int l_cdate = 1;
        const int l_rdate = 2;

        // Count the number of rows that matched per o_priority
        std::map<Data, Data> priority_count_real;
        for (int i = 0; i < rows_1; ++i) {
            if (orders_data[i][o_date] >= 30 && orders_data[i][o_date] < 70){
                for(int j = 0; j < rows_2; ++j){
                    if(orders_data[i][o_key] == line_data[j][l_key]
                        && line_data[j][l_cdate] < line_data[j][l_rdate]){

                        priority_count_real[orders_data[i][o_priority]]++;
                        break;

                    }
                }
            }
        }

        // Collect the result from secure execution
        std::map<Data, Data> priority_count_res;
        int rows_3 = output.getSize();
        for(int i = 0; i < rows_3; ++i){
            if(res_open_ptr[o_priority][i] != -1){
                priority_count_res[res_open_ptr[o_priority][i]] = res_open_ptr[o_sel][i];
            }
        }

        printf("Output count %ld==%ld\n", priority_count_real.size(), priority_count_res.size());

        return priority_count_real == priority_count_res;
    } else {
        return true;
    }
}


void planner_q4_test() {

    // Orders Table
    int orders_rows = 128;
    std::vector<std::string> orders_schema = {"[order_key]",
                                              "[order_date]", "[order_priority]",
                                              "[order_date-D1]", "[order_date-D2]",
                                              "[SEL]", "SEL"};
    Data **orders_data = getRandomData(orders_rows, orders_schema.size());
    for(int i = 0; i < orders_rows; ++i){
        // D1 = 30      D2 = 70
        orders_data[i][3] = orders_data[i][1] - 30;
        orders_data[i][4] = orders_data[i][1] - 70;
    }
    Relation orders = getTestRelation("ORDERS", orders_rows, orders_schema, orders_data);


    // Line table
    int line_row = 128;
    std::vector<std::string> line_schema = {"[order_key]",
                                            "[commit_date]", "[receipt_date]",
                                            "[commit_date-receipt_date]",
                                            "[SEL]"};
    Data **line_data = getRandomData(line_row, line_schema.size());
    for(int i = 0; i < line_row; ++i){
        line_data[i][3] = line_data[i][1] - line_data[i][2];
    }
    Relation line = getTestRelation("LINE", line_row, line_schema, line_data);


    //exchange seeds
    exchange_rsz_seeds();

    // QUERY OPTIMIZATION:
    orders["[SEL]"] = (orders["[order_date-D1]"] >= 0 ) & (orders["[order_date-D2]"] < 0);
    orders.bSelectionUsed = true;

    line["[SEL]"] = line["[commit_date]"] < line["[receipt_date]"];
    line.bSelectionUsed = true;

    ScanOperator orders_Scan(orders, 4);
    ScanOperator line_Scan(line, 4);
    Expression join_expr(ExpressionType::Equal,
                         Expression("ORDERS", "[order_key]"),
                         Expression("LINE", "[order_key]"));
    SemiJoinOperator semiJoin(orders_Scan, line_Scan, Expression("ORDERS", "[order_key]"), Expression("LINE", "[order_key]"), 2);
    Expression star(ExpressionType::Star);
    GroupByOperator groupBy(semiJoin, GroupByType::COUNT, {Expression("ORDERS", "[SEL]"), Expression("ORDERS", "[order_priority]")}, {Expression("ORDERS", "SEL"), Expression("ORDERS", "[SEL]")}, {GroupByType::COUNT}, {star}, true, 2);
    // NO ORDER BY because GROUP BY will sort while grouping

    // Execution
    groupBy.execute();

    bool test_res = checkQ4Correctness(orders_data, orders_rows, orders_schema,
                                        line_data, line_row, line_schema,
                                        groupBy.output);

    assert(test_res == true);

    // free temp data
    free(orders_data);
    free(line_data);
    orders.freeRelation();
    line.freeRelation();
}


int main(int argc, char **argv) {

    // initialize communication
    init(argc, argv);
    init_sharing();

    // Testing and timing the query
    timeTest(planner_q4_test, "planner_q4_test");

    // tear down communication
    MPI_Finalize();
    return 0;
}