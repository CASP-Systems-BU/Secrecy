#ifndef SECRECY_COLUMN_H
#define SECRECY_COLUMN_H

#include<iostream>

class Column {
public:
    virtual Column &operator+(const Column &other) const { return *new Column(); }
    virtual Column &operator-(const Column &other) const { return *new Column(); }
    virtual Column &operator*(const Column &other) const { return *new Column(); }
    virtual Column &operator-() const { return *new Column(); }

    virtual Column &operator+(const Data &other) const { return *new Column(); }
    virtual Column &operator-(const Data &other) const { return *new Column(); }
    virtual Column &operator*(const Data &other) const { return *new Column(); }


    virtual Column &operator^(const Column &other) const { return *new Column(); }
    virtual Column &operator&(const Column &other) const { return *new Column(); }
    virtual Column &operator|(const Column &other) const { return *new Column(); }
    virtual Column &operator~() const { return *new Column(); }
    virtual Column &operator!() const { return *new Column(); }

    virtual Column &operator^(const Data &other) const { return *new Column(); }
    virtual Column &operator&(const Data &other) const { return *new Column(); }
    virtual Column &operator|(const Data &other) const { return *new Column(); }

    virtual Column &operator==(const Column &other) const { return *new Column(); }
    virtual Column &operator!=(const Column &other) const { return *new Column(); }
    virtual Column &operator>(const Column &other) const { return *new Column(); }
    virtual Column &operator>=(const Column &other) const { return *new Column(); }
    virtual Column &operator<(const Column &other) const { return *new Column(); }
    virtual Column &operator<=(const Column &other) const { return *new Column(); }

    virtual Column &operator==(const Data &other) const { return *new Column(); }
    virtual Column &operator!=(const Data &other) const { return *new Column(); }
    virtual Column &operator>(const Data &other) const { return *new Column(); }
    virtual Column &operator>=(const Data &other) const { return *new Column(); }
    virtual Column &operator<(const Data &other) const { return *new Column(); }
    virtual Column &operator<=(const Data &other) const { return *new Column(); }

    virtual Column &operator=(const Column &other) { return *this; }

    virtual int getIndex() const {
        printf(" from base class");
        return -1;
    }

    virtual long long getTableAddress() const {
        return -1;
    }

    virtual Column &setTableAddress(long long address) {
        return *this;
    }
};

#endif //SECRECY_COLUMN_H
