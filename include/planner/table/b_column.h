#ifndef SECRECY_B_COLUMN_H
#define SECRECY_B_COLUMN_H

#include "derived_column.h"

// keep track of the initiated new AColumn
template<typename T, typename T2>
class BColumn : public DerivedColumn<T, T2> {

public:
    BColumn(Column &other) : DerivedColumn<T, T2>(other) {}
    BColumn(DerivedColumn<T, T2> &other) : DerivedColumn<T, T2>(other) {}

    // Comparison with other column
    Column &operator==(const Column &other) const {
        auto res = this->operatorFunction(select_EQ_array_b, *this, other);
        return *new BColumn<T, T2>(res);
    }

    Column &operator!=(const Column &other) const {
        Column &c = ~((*this) == other);
        return c;
    }

    Column &operator>(const Column &other) const {
        auto res = this->operatorFunction(select_GT_array_b, *this, other);
        return *new BColumn<T, T2>(res);
    }

    Column &operator>=(const Column &other) const {
        auto res = this->operatorFunction(select_GEQ_array_b, *this, other);
        return *new BColumn<T, T2>(res);
    }

    Column &operator<(const Column &other) const {
        auto res = this->operatorFunction(select_GT_array_b, other, *this);
        return *new BColumn<T, T2>(res);
    }

    Column &operator<=(const Column &other) const {
        auto res = this->operatorFunction(select_GEQ_array_b, other, *this);
        return *new BColumn<T, T2>(res);
    }

    // Boolean Operations
    Column &operator&(const Column &other) const {
        auto res = this->operatorFunction(select_AND_array_b, *this, other);
        return *new BColumn<T, T2>(res);
    }

    Column &operator|(const Column &other) const {
        auto res = this->operatorFunction(select_OR_array_b, *this, other);
        return *new BColumn<T, T2>(res);
    }

    Column &operator^(const Column &other) const {
        int rows = this->getSize();

        auto second_derived = reinterpret_cast<const DerivedColumn<T, T2> *>(&other);

        auto shares_1 = this->getShares();
        auto shares_2 = second_derived->getShares();

        auto shares_1_ptr = shares_1.get();
        auto shares_2_ptr = shares_2.get();

        DerivedColumn<T, T2> res_data(rows);
        auto res = res_data.getShares();
        auto res_ptr = res.get();

        for(int i = 0; i < rows; ++i){
            res_ptr[0][i] = shares_1_ptr[0][i] ^ shares_2_ptr[0][i];
            res_ptr[1][i] = shares_1_ptr[1][i] ^ shares_2_ptr[1][i];
        }

        return *new BColumn<T, T2>(res_data);
    }

    Column &operator~() const {
        auto res = this->operatorFunction(select_NOT_array_b, *this, *this);
        return *new BColumn<T, T2>(res);
    }

    Column &operator!() const {
        DerivedColumn<T, T2> res_data(this->getSize());
        auto res = BColumn(res_data);
        auto res_array = res.getShares();
        auto res_array_ptr = res_array.get();

        auto this_array = this->getShares();
        auto this_array_ptr = this_array.get();

        for(int i = 0; i < this->getSize(); ++i){
            res_array_ptr[0][i] = this_array_ptr[0][i] ^ 1;
            res_array_ptr[1][i] = this_array_ptr[1][i] ^ 1;
        }

        return *new BColumn<T, T2>(res);
    }

    Column &operator < (const Data &other) const {
        int rows = this->getSize();

        if(other != 0){
            DerivedColumn<T, T2> other_data(rows);
            BColumn other_array = BColumn(other_data);
            auto other_array_ptr  = other_array.getShares().get();
            if (get_rank() == 0){
                for(int i = 0; i < rows; ++i){
                    other_array_ptr[0][i] = other;
                    other_array_ptr[1][i] = 0;
                }
            } else if (get_rank() == 1){
                for(int i = 0; i < rows; ++i){
                    other_array_ptr[0][i] = 0;
                    other_array_ptr[1][i] = 0;
                }
            } else if (get_rank() == 2){
                for(int i = 0; i < rows; ++i){
                    other_array_ptr[0][i] = 0;
                    other_array_ptr[1][i] = other;
                }
            }
            return (*this) < other_array;
        }

        assert (other  == 0);

        auto shares_1 = this->getShares();
        auto shares_1_ptr = shares_1.get();

        auto res_shares = allocateColumnWiseContent<T>(rows, 1);
        auto res_shares_ptr = res_shares.get();

        // if the public constant is zero
        if (other == 0){
            for (int i = 0; i < rows; ++i) {
                res_shares_ptr[0][i] = ltz_b(shares_1_ptr[0][i]);
                res_shares_ptr[1][i] = ltz_b(shares_1_ptr[1][i]);
            }
        }

        DerivedColumn<T, T2> res (res_shares, rows, ColumnType::INTERMEDIATE_ARRAY);

        return *new BColumn<T, T2>(res);
    }

    Column &operator > (const Data &other) const {
        int rows = this->getSize();

        DerivedColumn<T, T2> other_data(rows);
        BColumn other_array = BColumn(other_data);
        auto other_array_ptr  = other_array.getShares().get();
        if (get_rank() == 0){
            for(int i = 0; i < rows; ++i){
                other_array_ptr[0][i] = other;
                other_array_ptr[1][i] = 0;
            }
        } else if (get_rank() == 1){
            for(int i = 0; i < rows; ++i){
                other_array_ptr[0][i] = 0;
                other_array_ptr[1][i] = 0;
            }
        } else if (get_rank() == 2){
            for(int i = 0; i < rows; ++i){
                other_array_ptr[0][i] = 0;
                other_array_ptr[1][i] = other;
            }
        }

        return (*this) > other_array;
    }

    Column &operator >= (const Data &other) const {
        return ~((*this) < other);
    }
};


#endif //SECRECY_B_COLUMN_H
