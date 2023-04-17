#ifndef SECRECY_COLUMN_CONVERSION_H
#define SECRECY_COLUMN_CONVERSION_H

#include "a_column.h"
#include "b_column.h"
#include "../../core/primitives.h"
#include "../../common/communication.h"

template<typename T, typename T2>
AColumn<T, T2> convertToAColumn(const Column &other) {
    auto other_ = reinterpret_cast<const DerivedColumn<BShare, BShareTable> *>(&other);
    auto other_shares = other_->getShares();
    auto other_shares_ptr = other_shares.get();

    int rows = other_->getSize();

    auto res_shares = allocateColumnWiseContent<T>(rows, 1);
    auto res_shares_ptr = res_shares.get();

    AShare *r_a = new AShare[rows];
    BShare *r_b = new BShare[rows];
    for (int i = 0; i < rows; ++i) {
        r_a[i] = get_next_r();
        r_b[i] = get_next_rb();
    }

    convert_single_bit_array(other_shares_ptr[0], r_a, r_b, rows, res_shares_ptr[0]);

    free(r_a);
    free(r_b);

    // communication round
    exchangeArrayWithParties<AShare>(res_shares_ptr[0], res_shares_ptr[1], rows, 1, PRIMITIVES_SHARE_TAG, 1);
    DerivedColumn<T, T2> res(res_shares, rows);

    return AColumn<T, T2>(res);
}


template<typename T, typename T2>
BColumn<T, T2> convertToBColumn(const Column &other) {
    auto other_ = reinterpret_cast<const DerivedColumn<AShare, AShareTable> *>(&other);
    auto other_shares = other_->getShares();
    auto other_shares_ptr = other_shares.get();

    int rows = other_->getSize();

    auto res_shares = allocateColumnWiseContent<T>(rows, 1);
    auto res_shares_ptr = res_shares.get();

    convert_a_to_b_array(other_shares_ptr[0], other_shares_ptr[1], res_shares_ptr[0], res_shares_ptr[1], rows);
    DerivedColumn<T, T2> res(res_shares, rows);

    return BColumn<T, T2>(res);
}


#endif //SECRECY_COLUMN_CONVERSION_H