#ifndef SECRECY_DERIVED_COLUMN_H
#define SECRECY_DERIVED_COLUMN_H

#include "column.h"
#include "../../common/table.h"
#include "../../core/wrappers/select.h"
#include "../../common/constants.h"

enum ColumnType {
    TABLE_REFERENCE,
    ARRAY_REFERENCE,
    INTERMEDIATE_ARRAY,
    NONE,
    DELETED
};

template<typename T, typename T2>
class DerivedColumn : public Column {
protected:
    ColumnType columnType;

    T2 *shareTable;
    int index;
    T *rowShares;
    int step;
    bool allocated = false;

    std::shared_ptr<T *> columnShares;
    int size;

    // An ith share for a column.
    DerivedColumn(T2 *shareTable, int index, T *rowShares, int step) :
            shareTable(shareTable), index(index), rowShares(rowShares),
            step(step) {
        size = 0;
        columnType = ColumnType::TABLE_REFERENCE;
        columnShares = nullptr;

    }

    // Copy constructor from Column
    DerivedColumn(Column &other) {
        auto other_ = reinterpret_cast<DerivedColumn<T, T2> *>(&other);

        // starting copying now
        this->columnType = other_->columnType;

        this->shareTable = other_->shareTable;
        this->index = other_->index;
        this->rowShares = other_->rowShares;
        this->step = other_->step;

        this->columnShares = other_->columnShares;
        this->size = other_->size;
    }

public:
    // TABLE_REFERENCE
    DerivedColumn(T2 *shareTable, int index) :
            shareTable(shareTable), index(index) {
        if (shareTable != nullptr) {
            columnType = ColumnType::TABLE_REFERENCE;
            size = shareTable->numRows;

            step = shareTable->numCols * CR_P.PARTY_REPLICATION;
        } else {
            columnType = ColumnType::NONE;
        }

        columnShares = nullptr;
    }

    void allocateRowShares(T2 *shareTable_ = nullptr, bool forceAllocation = false) {
        if (!allocated || forceAllocation) {
            if (columnType == ColumnType::TABLE_REFERENCE) {
                if (shareTable_ != nullptr) {
                    this->shareTable = shareTable_;
                }
                assert(shareTable->content != nullptr);
                int base_index = index * CR_P.PARTY_REPLICATION;
                rowShares = &shareTable->content[0][base_index];
            }
        }
        allocated = true;
    }

    // ARRAY_REFERENCE - INTERMEDIATE_ARRAY
    DerivedColumn(std::shared_ptr<T *> &columnShares, int size, ColumnType columnType = ARRAY_REFERENCE) :
            columnShares(columnShares), size(size), columnType(columnType) {

        rowShares = columnShares.get()[0];
        step = 1;
        index = -1;
    }

    DerivedColumn(int size, ColumnType columnType = ARRAY_REFERENCE) :
            size(size), columnType(columnType) {
        columnShares = allocateColumnWiseContent<T>(size, 1);
        auto columnSharesPtr = columnShares.get();

        rowShares = columnSharesPtr[0];
        step = 1;
        index = -1;
    }

    // Not valid DerivedColumn instance.
    DerivedColumn() : columnType(ColumnType::NONE) {}

    virtual int getIndex() const {
        if (columnType == ColumnType::TABLE_REFERENCE) {
            return index;
        } else {
            return -1;
        }
    }

    inline int getSize() const {
        if (columnType == ColumnType::TABLE_REFERENCE) {
            return shareTable->numRows;
        } else {
            return size;
        }
    }

    virtual long long getTableAddress() const {
        return (long long) shareTable;
    }

    virtual Column &setTableAddress(long long address) {
        shareTable = (T2 *) address;
        return *this;
    }

    inline T &operator[](const int &row) {
        return rowShares[row * step];
    }

    DerivedColumn getShare(int i = 0) const {
        if (columnType == ColumnType::TABLE_REFERENCE) {
            return DerivedColumn(shareTable, index, &rowShares[i], step);
        } else if (columnType == ColumnType::ARRAY_REFERENCE
                   || columnType == ColumnType::INTERMEDIATE_ARRAY) {
            return DerivedColumn(columnShares + i, size, columnType);
        } else {
            return DerivedColumn();
        }
    }

    std::shared_ptr<T *> getShares() const {
        switch (columnType) {
            case ColumnType::TABLE_REFERENCE:
                return getFromTable<T, T2>(*shareTable, index);
                break;
            case ColumnType::ARRAY_REFERENCE:
            case ColumnType::INTERMEDIATE_ARRAY:
                return columnShares;
                break;
            default:
                return nullptr;
                break;
        }
    }

    ColumnType getType() const {
        return columnType;
    }

    DerivedColumn<T, T2>
    operatorFunction(std::shared_ptr<T *>(func)(std::shared_ptr<T *>, std::shared_ptr<T *>, const int &),
                     const Column &first, const Column &other) const {
        auto first_derived = reinterpret_cast<const DerivedColumn *>(&first);
        auto second_derived = reinterpret_cast<const DerivedColumn *>(&other);

        auto shares_1 = first_derived->getShares();
        auto shares_2 = second_derived->getShares();

        auto shares_3 = func(shares_1, shares_2, size);

        return DerivedColumn<T, T2>(shares_3, size, ColumnType::INTERMEDIATE_ARRAY);
    }

    DerivedColumn &operator=(const Column &other) {
        auto other_ = reinterpret_cast<const DerivedColumn *>(&other);
        auto other__ = other_->getShares();
        (*this) = other__;
        return (*this);
    }

    DerivedColumn &operator=(const DerivedColumn &other) {
        const std::shared_ptr<T *> other_ = other.getShares();
        (*this) = other_;
        return (*this);
    }

    DerivedColumn &operator=(const std::shared_ptr<T *> &other) {
        auto otherArray_ptr = other.get();

        int col;
        switch (columnType) {
            case ColumnType::TABLE_REFERENCE:
                allocateRowShares();
                col = index * CR_P.PARTY_REPLICATION;
                for (int i = 0; i < size; ++i) {
                    for (int j = 0; j < CR_P.PARTY_REPLICATION; ++j) {
                        shareTable->content[i][col + j] = otherArray_ptr[j][i];
                    }
                }
                break;
            case ColumnType::ARRAY_REFERENCE: {
                auto this_ptr = columnShares.get();
                for (int j = 0; j < CR_P.PARTY_REPLICATION; ++j) {
                    for (int i = 0; i < this->size; ++i) {
                        this_ptr[j][i] = otherArray_ptr[j][i];
                    }
                }
            }
                break;
            case ColumnType::INTERMEDIATE_ARRAY:
                columnShares = other;
                break;
            default:
                break;
        }
        return *this;
    }
};

#endif //SECRECY_DERIVED_COLUMN_H
