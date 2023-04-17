#include "../../../include/core/wrappers/select.h"
#include "../../../include/common/communication.h"
#include "../../../include/common/table.h"

// Assumes a computation without communication layer for each operator
// TODO: inline use functions and make input arguments const
// TODO: create function for creating (local - remote) rowShares
// TODO: create an additional optional parameter to reuse allocated memory

//////////////////////////////////////////////////////////
///       Comparison        /      (Boolean - table)     /
//////////////////////////////////////////////////////////
// TODO: use native primitives instead of reusing the arrays implementation
void select_EQ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo,
                       BShareTable &output, const int &colOne,
                       const int &colTwo, const int &colOut) {

    // Input validation
    validateTableInput<BShareTable>(inputOne, inputTwo, output, colOne, colTwo, colOut);

    // Reduce to array primitives
    auto shares_1 = getFromTable<BShare, BShareTable>(inputOne, colOne);
    auto shares_2 = getFromTable<BShare, BShareTable>(inputTwo, colTwo);

    // Computation - no for Communication round because it is done inside this call
    auto shares = select_EQ_array_b(shares_1, shares_2, inputOne.numRows);

    addToTable<BShare, BShareTable>(output, shares, colOut);
}

// TODO: use native primitives instead of reusing the arrays implementation
void select_GT_table_b(const BShareTable &inputOne, const BShareTable &inputTwo,
                       BShareTable &output, const int &colOne, const int &colTwo,
                       const int &colOut) {

    // Input validation
    validateTableInput<BShareTable>(inputOne, inputTwo, output, colOne, colTwo, colOut);

    // Reduce to array primitives
    auto shares_1 = getFromTable<BShare, BShareTable>(inputOne, colOne);
    auto shares_2 = getFromTable<BShare, BShareTable>(inputTwo, colTwo);

    // Computation - no for Communication round because it is done inside this call
    auto shares = select_GT_array_b(shares_1, shares_2, inputOne.numRows);

    addToTable<BShare, BShareTable>(output, shares, colOut);
}

// TODO: use native primitives instead of reusing the arrays implementation
void select_GEQ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo,
                        BShareTable &output, const int &colOne, const int &colTwo,
                        const int &colOut) {

    // Input validation
    validateTableInput<BShareTable>(inputOne, inputTwo, output, colOne, colTwo, colOut);

    // Reduce to array primitives
    auto shares_1 = getFromTable<BShare, BShareTable>(inputOne, colOne);
    auto shares_2 = getFromTable<BShare, BShareTable>(inputTwo, colTwo);

    // Computation - no for Communication round because it is done inside this call
    auto shares = select_GEQ_array_b(shares_1, shares_2, inputOne.numRows);

    addToTable<BShare, BShareTable>(output, shares, colOut);
}

// TODO: is this the equality for pre public constant?
// TODO: needs a preprocessing step (c-att) & (att-c)
void select_EQZ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo,
                        BShareTable &output, const int &colOne,
                        const int &colTwo, const int &colOut) {

    // Input validation
    validateTableInput<BShareTable>(inputOne, inputTwo, output, colOne, colTwo, colOut);

    // Output is temporarily in an array for easy sharing communication
    auto shares = allocateColumnWiseContent<BShare>(output.numRows, 1);
    auto shares_ptr = shares.get()[0];

    // Computation
    // r = (c-att < 0) ^ (att-c < 0) ^ 1
    BShare z1, z2;
    for (int i = 0; i < inputOne.numRows; i++) {
        z1 = ltz_b(inputOne.content[i][colOne]);
        z2 = ltz_b(inputTwo.content[i][colTwo]);
        shares_ptr[i] = z1 ^ z2 ^ 1;
    }

    // Communication round - add rowShares to table - free temp arrays
    exchangeAndAddToTable(output, shares, colOut);
}

// TODO: is this the greater than for pre public constant?
// TODO: needs a preprocessing step (c-att)
void select_GTZ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo,
                        BShareTable &output, const int &colOne, const int &colTwo,
                        const int &colOut) {

    // Input validation
    validateTableInput<BShareTable>(inputOne, output, colOne, colOut);

    // Output is temporarily in an array for easy sharing communication
    auto shares = allocateColumnWiseContent<BShare>(output.numRows, 1);
    auto shares_ptr = shares.get()[0];

    // Computation
    // r = (c-att < 0)
    for (int i = 0; i < inputOne.numRows; i++) {
        shares_ptr[i] = ltz_b(inputOne.content[i][colOne]);
    }

    // Communication round - add rowShares to table - free temp arrays
    exchangeAndAddToTable(output, shares.get(), colOut);
}

// TODO: is this the greater than or equal for pre public constant?
// TODO: needs a preprocessing step (att - c)
void select_GEQZ_table_b(const BShareTable &inputOne, const BShareTable &inputTwo,
                         BShareTable &output, const int &colOne, const int &colTwo,
                         const int &colOut) {

    // Input validation
    validateTableInput<BShareTable>(inputOne, output, colOne, colOut);

    // Output is temporarily in an array for easy sharing communication
    auto shares = allocateColumnWiseContent<BShare>(output.numRows, 1);
    auto shares_ptr = shares.get()[0];

    // Computation
    // r = ~((att - c) < 0)
    for (int i = 0; i < inputOne.numRows; i++) {
        shares_ptr[i] = (ltz_b(inputOne.content[i][colOne]) & 1) ^ 1;
    }

    // Communication round - add rowShares to table - free temp arrays
    exchangeAndAddToTable(output, shares.get(), colOut);
}

//////////////////////////////////////////////////////////
///       Comparison        /      (Boolean - array)     /
////////////////////////////////////////////////////////// 

ColumnData select_EQ_array_b(ColumnData inputOne, ColumnData inputTwo,
                             const int &rows) {

    // allocate output
    auto shares = allocateColumnWiseContent<BShare>(rows, 1);

    // Computation
    eq_b_array(inputOne.get()[0], inputOne.get()[1], inputTwo.get()[0], inputTwo.get()[1], rows, shares.get()[0]);

    // Communication round
    exchangeArrayWithParties<BShare>(shares.get()[0], shares.get()[1], rows, 1, PRIMITIVES_SHARE_TAG, 1);

    return shares;
}

ColumnData select_GT_array_b(ColumnData inputOne, ColumnData inputTwo,
                             const int &rows) {

    // allocate output
    auto shares = allocateColumnWiseContent<BShare>(rows, 1);
    auto shares_ptr = shares.get()[0];
    auto shares_bit = allocateColumnWiseContent<BitShare>(rows, 1);
    auto shares_bit_ptr = shares_bit.get()[0];


    // Computation
    // TODO: instead of casting, unify the format for all similar functions
    greater_batch(inputOne.get()[0], inputOne.get()[1], inputTwo.get()[0], inputTwo.get()[1], rows,
                  shares_bit.get()[0]);

    for (int i = 0; i < rows; ++i) {
        shares_ptr[i] = shares_bit_ptr[i] & 1;
    }

    // Communication round
    exchangeArrayWithParties<BShare>(shares.get()[0], shares.get()[1], rows, 1, PRIMITIVES_SHARE_TAG, 1);

    return shares;
}

ColumnData select_GEQ_array_b(ColumnData inputOne, ColumnData inputTwo,
                              const int &rows) {

    // allocate output
    auto shares = allocateColumnWiseContent<BShare>(rows, 1);
    auto shares_ptr = shares.get()[0];
    auto shares_bit = allocateColumnWiseContent<BitShare>(rows, 1);
    auto shares_bit_ptr = shares_bit.get()[0];

    // Computation
    // TODO: instead of casting, unify the format for all similar functions
    geq_batch(inputOne.get()[0], inputOne.get()[1], inputTwo.get()[0], inputTwo.get()[1], rows, shares_bit.get()[0]);

    for (int i = 0; i < rows; ++i) {
        shares_ptr[i] = shares_bit_ptr[i] & 1;
    }

    // Communication round
    exchangeArrayWithParties<BShare>(shares.get()[0], shares.get()[1], rows, 1, PRIMITIVES_SHARE_TAG, 1);

    return shares;
}


ColumnData select_EQZ_array_b(ColumnData inputOne, ColumnData inputTwo,
                              const int &rows) {

    // allocate output
    auto shares = allocateColumnWiseContent<BShare>(rows, 1);
    auto shares_ptr = shares.get()[0];
    auto inputOne_ptr = inputOne.get()[0];
    auto inputTwo_ptr = inputTwo.get()[0];

    // Computation
    // r = (c-att < 0) ^ (att-c < 0) ^ 1
    BShare z1, z2;
    for (int i = 0; i < rows; i++) {
        z1 = ltz_b(inputOne_ptr[i]);
        z2 = ltz_b(inputTwo_ptr[i]);
        shares_ptr[i] = z1 ^ z2 ^ 1;
    }

    // Communication round
    exchangeArrayWithParties<BShare>(shares.get()[0], shares.get()[1], rows, 1, PRIMITIVES_SHARE_TAG, 1);

    return shares;
}

ColumnData select_GTZ_array_b(ColumnData inputOne, ColumnData inputTwo,
                              const int &rows) {

    // allocate output
    auto shares = allocateColumnWiseContent<BShare>(rows, 1);
    auto shares_ptr = shares.get()[0];
    auto inputOne_ptr = inputOne.get()[0];


    // Computation
    // r = (c-att < 0)
    for (int i = 0; i < rows; i++) {
        shares_ptr[i] = ltz_b(inputOne_ptr[i]);
    }

    // Communication round
    exchangeArrayWithParties<BShare>(shares.get()[0], shares.get()[1], rows, 1, PRIMITIVES_SHARE_TAG, 1);

    return shares;
}

ColumnData select_GEQZ_array_b(ColumnData inputOne, ColumnData inputTwo,
                               const int &rows) {

    // allocate output
    auto shares = allocateColumnWiseContent<BShare>(rows, 1);
    auto shares_ptr = shares.get()[0];
    auto inputOne_ptr = inputOne.get()[0];

    // Computation
    // r = ~((att - c) < 0)
    for (int i = 0; i < rows; i++) {
        shares_ptr[i] = (ltz_b(inputOne_ptr[i]) & 1) ^ 1;
    }

    // Communication round
    exchangeArrayWithParties<BShare>(shares.get()[0], shares.get()[1], rows, 1, PRIMITIVES_SHARE_TAG, 1);

    return shares;
}


//////////////////////////////////////////////////////////
///   Boolean Operations    /      (Boolean - table)     /
//////////////////////////////////////////////////////////
void select_AND_table_b(const BShareTable &inputOne, const BShareTable &inputTwo,
                        BShareTable &output, const int &colOne,
                        const int &colTwo, const int &colOut) {

    // Input validation
    validateTableInput<BShareTable>(inputOne, inputTwo, output, colOne, colTwo, colOut);

    int rows = inputOne.numRows;
    BShares randomNumbers = new BShare[rows];
    get_next_rb_array(randomNumbers, rows);


    // Output is temporarily in an array for easy sharing communication
    auto shares = allocateColumnWiseContent<BShare>(output.numRows, 1);

    auto shares_ptr = shares.get();
    auto shares_ptr_ = shares_ptr[0];

    // Computation
    for (int i = 0; i < inputOne.numRows; i++) {
        shares_ptr_[i] = and_b(inputOne.content[i][colOne],
                               inputOne.content[i][colOne + 1],
                               inputTwo.content[i][colTwo],
                               inputTwo.content[i][colTwo + 1], randomNumbers[i]);
    }

    // Free temp array
    free(randomNumbers);

    // Communication round - add rowShares to table - free temp arrays
    exchangeAndAddToTable(output, shares_ptr, colOut);
}


//////////////////////////////////////////////////////////
///   Boolean Operations    /      (Boolean - array)     /
////////////////////////////////////////////////////////// 
ColumnData select_AND_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows) {
    // first generate the random numbers
    BShares randomNumbers = new BShare[rows];
    get_next_rb_array(randomNumbers, rows);

    // allocate memory for result;
    auto shares = allocateColumnWiseContent<BShare>(rows, 1);

    // Start MPC protocol
    and_b_array(inputOne.get()[0], inputOne.get()[1], inputTwo.get()[0], inputTwo.get()[1], randomNumbers, rows,
                shares.get()[0]);

    // Free temp arrays
    free(randomNumbers);

    // now a round of communication
    exchangeArrayWithParties<BShare>(shares.get()[0], shares.get()[1], rows, 1, PRIMITIVES_SHARE_TAG, 1);

    return shares;
}

// Using and and de morgan rule.
ColumnData select_OR_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows) {
    // First negate the inputs
    auto inputOne_ = select_NOT_array_b(inputOne, nullptr, rows);
    auto inputTwo_ = select_NOT_array_b(inputTwo, nullptr, rows);

    // Use the select_AND_array_b - Free temp inputs
    auto output_ = select_AND_array_b(inputOne_, inputTwo_, rows);

    // Negate the output - Free temp output
    auto output = select_NOT_array_b(output_, nullptr, rows);

    return output;
}

// TODO: just negating one of the replicated rowShares, right?
ColumnData select_NOT_array_b(ColumnData inputOne, ColumnData inputTwo, const int &rows) {

    // allocate memory for result;
    auto shares = allocateColumnWiseContent<BShare>(rows, 1);
    auto shares_ptr = shares.get();
    auto input_ptr = inputOne.get();

    // Start MPC protocol
    if (get_rank() == 0) {
        for (int i = 0; i < rows; ++i) {
            shares_ptr[0][i] = input_ptr[0][i];
            shares_ptr[1][i] = ~input_ptr[1][i];
        }
    } else if (get_rank() == 1) {
        for (int i = 0; i < rows; ++i) {
            shares_ptr[0][i] = ~input_ptr[0][i];
            shares_ptr[1][i] = input_ptr[1][i];
        }
    } else {
        for (int i = 0; i < rows; ++i) {
            shares_ptr[0][i] = input_ptr[0][i];
            shares_ptr[1][i] = input_ptr[1][i];
        }
    }

    // now a round of communication
    // No communication for the not

    return shares;
}

//////////////////////////////////////////////////////////
///  Arithmetic Operations  /    (Arithmetic - table)    /
////////////////////////////////////////////////////////// 
void select_ADD_table(const AShareTable &inputOne, const AShareTable &inputTwo,
                      AShareTable &output, const int &colOne,
                      const int &colTwo, const int &colOut) {

    // Input validation
    validateTableInput<AShareTable>(inputOne, inputTwo, output, colOne, colTwo, colOut);

    int col_11 = colOut;
    int col_12 = col_11 + 1;

    int col_21 = colOut;
    int col_22 = col_21 + 1;

    int col_31 = colOut;
    int col_32 = col_31 + 1;

    // Computation
    for (int i = 0; i < inputOne.numRows; ++i) {
        output.content[i][col_31] = inputOne.content[col_11][i] + inputTwo.content[col_21][i];
        output.content[i][col_32] = inputOne.content[col_12][i] + inputTwo.content[col_22][i];
    }

    // No communication round
}

void select_SUB_table(const AShareTable &inputOne, const AShareTable &inputTwo,
                      AShareTable &output, const int &colOne,
                      const int &colTwo, const int &colOut) {

    // Input validation
    validateTableInput<AShareTable>(inputOne, inputTwo, output, colOne, colTwo, colOut);

    int col_11 = colOut;
    int col_12 = col_11 + 1;

    int col_21 = colOut;
    int col_22 = col_21 + 1;

    int col_31 = colOut;
    int col_32 = col_31 + 1;

    // Computation
    for (int i = 0; i < inputOne.numRows; ++i) {
        output.content[i][col_31] = inputOne.content[col_11][i] - inputTwo.content[col_21][i];
        output.content[i][col_32] = inputOne.content[col_12][i] - inputTwo.content[col_22][i];
    }

    // No communication round
}

// TODO: batch the random Ashare of zero generation
void select_MULT_table(const AShareTable &inputOne, const AShareTable &inputTwo,
                       AShareTable &output, const int &colOne,
                       const int &colTwo, const int &colOut) {

    // Input validation
    validateTableInput<AShareTable>(inputOne, inputTwo, output, colOne, colTwo, colOut);

    // Output is temporarily in an array for easy sharing communication
    auto shares = allocateColumnWiseContent<AShare>(output.numRows, 1);
    auto shares_ptr = shares.get()[0];

    int col_11 = colOut;
    int col_12 = col_11 + 1;

    int col_21 = colOut;
    int col_22 = col_21 + 1;


    // Computation
    for (int i = 0; i < inputOne.numRows; i++) {
        shares_ptr[i] = mul(inputOne.content[col_11][i], inputOne.content[col_12][i],
                            inputTwo.content[col_21][i], inputTwo.content[col_21][i],
                            get_next_r());
    }


    // Communication round - add rowShares to table - free temp arrays
    exchangeAndAddToTable<AShare, AShareTable>(output, shares.get(), colOut);
}

//////////////////////////////////////////////////////////
///  Arithmetic Operations  /    (Arithmetic - array)    /
////////////////////////////////////////////////////////// 
ColumnData select_ADD_array_a(ColumnData inputOne, ColumnData inputTwo, const int &rows) {
    // Output
    auto shares = allocateColumnWiseContent<AShare>(rows, 1);
    auto shares_ptr = shares.get();
    auto input1_ptr = inputOne.get();
    auto input2_ptr = inputTwo.get();

    // Computation
    for (int i = 0; i < rows; ++i) {
        shares_ptr[0][i] = input1_ptr[0][i] + input2_ptr[0][i];
        shares_ptr[1][i] = input1_ptr[1][i] + input2_ptr[1][i];
    }

    // No communication round
    return shares;
}

ColumnData select_SUB_array_a(ColumnData inputOne, ColumnData inputTwo, const int &rows) {
    // Output 
    auto shares = allocateColumnWiseContent<AShare>(rows, 1);
    auto shares_ptr = shares.get();

    auto input1_ptr = inputOne.get();
    auto input2_ptr = inputTwo.get();

    // Computation
    for (int i = 0; i < rows; ++i) {
        shares_ptr[0][i] = input1_ptr[0][i] - input2_ptr[0][i];
        shares_ptr[1][i] = input1_ptr[1][i] - input2_ptr[1][i];
    }

    // No communication round
    return shares;
}

// TODO: batch the random Ashare of zero generation
ColumnData select_MULT_array_a(ColumnData inputOne, ColumnData inputTwo, const int &rows) {
    // Output 
    auto shares = allocateColumnWiseContent<AShare>(rows, 1);
    auto shares_ptr = shares.get();

    auto input1_ptr = inputOne.get();
    auto input2_ptr = inputTwo.get();


    // Computation
    for (int i = 0; i < rows; i++) {
        shares_ptr[0][i] = mul(input1_ptr[0][i], input1_ptr[1][i],
                               input2_ptr[0][i], input2_ptr[1][i],
                               get_next_r());
    }

    // now a round of communication
    exchangeArrayWithParties<AShare>(shares_ptr[0], shares_ptr[1], rows, 1, PRIMITIVES_SHARE_TAG, 1);

    return shares;
}