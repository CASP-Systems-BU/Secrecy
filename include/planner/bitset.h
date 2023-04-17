#ifndef SECRECY_BITSET_H
#define SECRECY_BITSET_H

// Contains information about relations included in problem solution
class BitSet {
public:
   // The type of the value
   typedef unsigned long value_t;
private:
   // The first bit
   static const value_t one = 1;

   // Constructor
   explicit BitSet(value_t mask) : mask(mask) { }

public:
    // The value
    value_t mask;
    // Constructor
    BitSet() : mask(0) { }
    // Sets a specific entry
    void set(unsigned i) { mask |= (one << i); }
    // Sets a specific entry
    void set(BitSet &m) { mask = m.mask; }
    // Checks if the given id is included in the bitset
    bool contains(unsigned i) const { return mask & (one << i); }

    // Equal operator
    bool operator==(const BitSet &o) const { return mask == o.mask; }
    // Overlaps with
    bool overlapsWith(const BitSet &o) const { return mask & o.mask; }
    // Union
    BitSet unionWith(const BitSet &o) const { return BitSet(mask | o.mask); }
};

#endif  // SECRECY_BITSET_H