#ifndef SECRECY_A_COLUMN_H
#define SECRECY_A_COLUMN_H

#include "derived_column.h"

// TODO: convert from BShare?
// keep track of the initiated new AColumn
template<typename T, typename T2>
class AColumn : public DerivedColumn<T, T2> {
public:
    // TODO: what happens on deleting the initial?
    //  >> check std::shared_ptr<T*> flow
    AColumn(Column &other) : DerivedColumn<T, T2>(other) {}
    AColumn(DerivedColumn<T, T2> &other) : DerivedColumn<T, T2>(other) {}

    // Arithmetic Operations
    Column &operator+(const Column &other) const {
        auto res = this->operatorFunction(select_ADD_array_a, *this, other);
        return *new AColumn<T, T2>(res);
    }

    Column &operator-(const Column &other) const {
        auto res = this->operatorFunction(select_SUB_array_a, *this, other);
        return *new AColumn<T, T2>(res);
    }

    Column &operator*(const Column &other) const {
        auto res = this->operatorFunction(select_MULT_array_a, *this, other);
        return *new AColumn<T, T2>(res);
    }

    Column & sum()const{
        std::shared_ptr<T *> res =  allocateColumnWiseContent<AShare>(1, 1);
        auto res_ptr = res.get();

        int rows = this->getSize();
        auto shares_1 = this->getShares();
        auto shares_1_ptr = shares_1.get();

        for (int i = 0; i < rows; ++i) {
            res_ptr[0][0] += shares_1_ptr[0][i];
            res_ptr[1][0] += shares_1_ptr[1][i];
        }

        auto res_intermediate = DerivedColumn<T, T2>(res, 1, ColumnType::INTERMEDIATE_ARRAY);

        return *new AColumn<T, T2>(res_intermediate);
    }

    // The public constant is secret shared as [other, 0, 0]
    Column &operator+(const Data &other) const {
        int rows = this->getSize();

        auto shares_1 = this->getShares();
        auto shares_1_ptr = shares_1.get();

        auto res_shares = allocateColumnWiseContent<AShare>(rows, 1);
        auto res_shares_ptr = res_shares.get();

        if (get_rank() == 0) {
            for (int i = 0; i < rows; ++i) {
                res_shares_ptr[0][i] = shares_1_ptr[i][0] + other;
                res_shares_ptr[1][i] = shares_1_ptr[i][1];
            }
        }else if(get_rank() == 1){
            for (int i = 0; i < rows; ++i) {
                res_shares_ptr[0][i] = shares_1_ptr[i][0];
                res_shares_ptr[1][i] = shares_1_ptr[i][1];
            }
        }else if(get_rank() == 2){
            for (int i = 0; i < rows; ++i) {
                res_shares_ptr[0][i] = shares_1_ptr[i][0];
                res_shares_ptr[1][i] = shares_1_ptr[i][1] + other;
            }
        }

        DerivedColumn<T, T2> res (res_shares, rows, ColumnType::INTERMEDIATE_ARRAY);

        return * new AColumn<T, T2>(res);
    }

    Column &operator-(const Data &other) const {
        return (*this) + (-other);
    }
};


#endif //SECRECY_A_COLUMN_H
