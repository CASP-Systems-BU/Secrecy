#include <stdio.h>
#include <assert.h>

#include "../test-utils.h"

void columnTest() {

    // // Table schema
    // int rows = 10000;
    // std::vector<std::string> r_a = {"x", "y", "z", "x+y", "x-y", "x*y", "x+y+z"};
    // std::vector<std::string> r_b = {"[x]", "[y]", "[z]", "[~x]", "[x&y]", "[x|y]", "[x==y]", "[x!=y]", "[x>y]", "[x<y]",
    //                                 "[x>=y]", "[x<=y]"};

    // // // Random data
    // Data **r_a_data = getRandomData(rows, r_a.size());
    // Data **r_b_data = getRandomData(rows, r_b.size());

    // // // Get rowShares for random data
    // Relation r = getTestRelation("R1", rows, r_a, r_a_data, r_b, r_b_data);

    // //exchange seeds
    // exchange_rsz_seeds();

    // // Arithmetic
    // // "x", "y", "z", "x+y", "x-y", "x*y", "x+y+z"
    // r["x+y"] = r["x"] + r["y"];
    // r["x-y"] = r["x"] - r["y"];
    // r["x*y"] = r["x"] * r["y"];
    // r["x+y+z"] = r["x"] + r["y"] + r["z"];


    // // Boolean
    // // "x", "y", "z", "~x" ,"x&y", "x|y", "x==y", "x!=y", ,"x>y", "x<y", "x>=y", "x<=y"}
    // r["[~x]"] = ~r["[x]"];
    // r["[x&y]"] = r["[x]"] & r["[y]"];
    // r["[x|y]"] = r["[x]"] | r["[y]"];
    // r["[x==y]"] = r["[x]"] == r["[y]"];
    // r["[x!=y]"] = r["[x]"] != r["[y]"];
    // r["[x>y]"] = r["[x]"] > r["[y]"];
    // r["[x<y]"] = r["[x]"] < r["[y]"];
    // r["[x>=y]"] = r["[x]"] >= r["[y]"];
    // r["[x<=y]"] = r["[x]"] <= r["[y]"];

    // // printf("%d res shares: %lld=%lld-%lld\t\t\t%lld=%lld-%lld\n", get_rank(), r.aShareTable.content[1][8], r.aShareTable.content[1][0], r.aShareTable.content[1][2], r.aShareTable.content[1][9], r.aShareTable.content[1][1], r.aShareTable.content[1][3]);
    // // printf("%d res shares: %lld=%lld>%lld\t\t\t%lld=%lld>%lld\n", get_rank(), r.bShareTable.content[0][16], r.bShareTable.content[0][0], r.bShareTable.content[0][2], r.bShareTable.content[0][17], r.bShareTable.content[0][1], r.bShareTable.content[0][3]);
    // // printf("%d res shares: %lld=~%lld\t\t\t%lld=~%lld\n", get_rank(), r.bShareTable.content[0][6], r.bShareTable.content[0][0], r.bShareTable.content[0][7], r.bShareTable.content[0][1]);


    // // Results assertion
    // // Column wise output
    // auto r_a_open_ = openTable(r.aShareTable);
    // auto r_a_open = r_a_open_.get();
    // auto r_b_open_ = openTable(r.bShareTable);
    // auto r_b_open = r_b_open_.get();

    // if (get_rank() == 0) {
    //     for (int i = 0; i < 10 && i < rows; ++i) {
    //         // checking inputs
    //         // printf("%lld==%lld\t\t%lld==%lld\t\t%lld==%lld\n", r_a_open[0][i], r_a_data[i][0], r_a_open[1][i], r_a_data[i][1], r_a_open[2][i], r_a_data[i][2]);
    //         assert(r_a_open[0][i] == r_a_data[i][0]);
    //         assert(r_a_open[1][i] == r_a_data[i][1]);
    //         assert(r_a_open[2][i] == r_a_data[i][2]);
    //         // printf("%lld==%lld\t\t%lld==%lld\t\t%lld==%lld\n", r_b_open[0][i], r_b_data[i][0], r_b_open[1][i], r_b_data[i][1], r_b_open[2][i], r_b_data[i][2]);
    //         assert(r_b_open[0][i] == r_b_data[i][0]);
    //         assert(r_b_open[1][i] == r_b_data[i][1]);
    //         assert(r_b_open[2][i] == r_b_data[i][2]);

    //         // printf("%lld=%lld+%lld\n", r_a_open[3][i], r_a_open[0][i], r_a_open[1][i]);
    //         assert(r_a_open[3][i] == r_a_open[0][i] + r_a_open[1][i]);

    //         // printf("%lld=%lld-%lld\n", r_a_open[4][i], r_a_open[0][i], r_a_open[1][i]);
    //         assert(r_a_open[4][i] == r_a_open[0][i] - r_a_open[1][i]);

    //         // printf("%lld=%lld*%lld\n", r_a_open[5][i], r_a_open[0][i], r_a_open[1][i]);
    //         assert(r_a_open[5][i] == r_a_open[0][i] * r_a_open[1][i]);

    //         // printf("%lld=%lld+%lld+%lld\n", r_a_open[6][i], r_a_open[0][i], r_a_open[1][i], r_a_open[2][i]);
    //         assert(r_a_open[6][i] == r_a_open[0][i] + r_a_open[1][i] + r_a_open[2][i]);

    //         // printf("%lld=~%lld\n", r_b_open[3][i], r_b_open[0][i]);
    //         assert(r_b_open[3][i] == (~r_b_open[0][i]));

    //         // printf("%lld=%lld&%lld\n", r_b_open[4][i], r_b_open[0][i], r_b_open[1][i]);
    //         assert(r_b_open[4][i] == (r_b_open[0][i] & r_b_open[1][i]));

    //         // printf("%lld=%lld|%lld\n", r_b_open[5][i], r_b_open[0][i], r_b_open[1][i]);
    //         assert(r_b_open[5][i] == (r_b_open[0][i] | r_b_open[1][i]));

    //         // printf("%lld=%lld==%lld\n", r_b_open[6][i], r_b_open[0][i], r_b_open[1][i]);
    //         assert((r_b_open[6][i] & 1) == (r_b_open[0][i] == r_b_open[1][i]));

    //         // TODO: mask the output with 1
    //         // printf("%lld=%lld!=%lld\n", (r_b_open[7][i]&1), r_b_open[0][i], r_b_open[1][i]);
    //         assert((r_b_open[7][i] & 1) == (r_b_open[0][i] != r_b_open[1][i]));

    //         // printf("%lld=%lld>%lld\t\t%d\n",r_b_open[8][i], r_b_open[0][i], r_b_open[1][i], r_b_open[0][i] > r_b_open[1][i]);
    //         assert(r_b_open[8][i] == (r_b_open[0][i] > r_b_open[1][i]));

    //         // printf("%lld=%lld<%lld\n", r_b_open[9][i], r_b_open[0][i], r_b_open[1][i]);
    //         assert(r_b_open[9][i] == (r_b_open[0][i] < r_b_open[1][i]));

    //         // printf("%lld=%lld>=%lld\n", r_b_open[10][i], r_b_open[0][i], r_b_open[1][i]);
    //         assert(r_b_open[10][i] == (r_b_open[0][i] >= r_b_open[1][i]));

    //         // printf("%lld=%lld<=%lld\n", r_b_open[11][i], r_b_open[0][i], r_b_open[1][i]);
    //         assert(r_b_open[11][i] == (r_b_open[0][i] <= r_b_open[1][i]));
    //     }
    // }

    // // free temp data
    // free(r_a_data);
    // free(r_b_data);
    // r.freeRelation();
}


int main(int argc, char **argv) {

    // initialize communication
    init(argc, argv);
    init_sharing();

    // Testing and timing the query
    timeTest(columnTest, "columnTest");

    // tear down communication
    MPI_Finalize();
    return 0;
}