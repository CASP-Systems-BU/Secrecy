#ifndef SECRECY_COSTMODEL_H
#define SECRECY_COSTMODEL_H

#include <cmath>
#include <iostream>

#define LENGTH 64           // Length of the share representation in bits
#define BATCH_SIZE 4096     // Max number of secret-shares exchanged between two parties per round
#define ALPHA 1             // Cost of local operation
#define BETA 8              // Cost ratio of remote to local operation

// The MPC cost model used by the plan generator
class CostModel {
public:
    // The cardinality type
    typedef double card_t;
    // Cost data type
    typedef double cost_t;
    // The actual MPC costs
    struct Cost {
        cost_t operation;
        cost_t synchronization;
        cost_t compOp;
        cost_t compSync;

        // Empty constructor
        Cost() : operation(0), synchronization(0), compOp(0), compSync(0) { };

        // Constructor
        Cost(cost_t operation, cost_t synchronization, cost_t compOp=0, cost_t compSync=0)
            : operation(operation), synchronization(synchronization), compOp(compOp), compSync(compSync) { };

        // Copy constructor
        Cost(const Cost &other) : operation(other.operation), synchronization(other.synchronization),
                                  compOp(other.compOp), compSync(other.compSync) { };

        // Addition in place
        Cost operator+=(const Cost &other) {
            this->operation += other.operation;
            this->synchronization += other.synchronization;
            this->compOp += other.compOp;
            this->compSync += other.compSync;
            return *this;
        }

        // Copy assignment
        Cost operator=(const Cost &other) {
            this->operation = other.operation;
            this->synchronization = other.synchronization;
            this->compOp = other.compOp;
            this->compSync = other.compSync;
            return *this;
        }

        // Addition
        Cost operator+(const Cost &other) const {
            Cost c = *this;
            c.operation += other.operation;
            c.synchronization += other.synchronization;
            c.compOp += other.compOp;
            c.compSync += other.compSync;
            return c;
        }

        // Less-than
        bool operator<(const Cost &other) const {
            if (this->operation < other.operation &&
                this->synchronization < other.synchronization &&
                this->compOp < other.compOp &&
                this->compSync < other.compSync)
                    return true;
            return false;
        }

        // Less-than-or-equal
        bool operator<=(const Cost &other) const {
            if (this->operation <= other.operation &&
                this->synchronization <= other.synchronization &&
                this->compOp <= other.compOp &&
                this->compSync <= other.compSync)
                return true;
            return false;
        }
    };

    // Oblivious AND cost
    static Cost andCost(bool bShare=true) {
        if (bShare) {
            Cost c = Cost(1,1);
            return c;
        }
        Cost c = a2bCost();
        c.operation++;
        c.synchronization++;
        return c;
    }

    // Oblivious XOR cost
    static Cost xorCost(bool bShare=true) {
        if (bShare) {
            Cost c = Cost(1,0);
            return c;
        }
        Cost c = a2bCost();
        c.operation++;
        return c;
    }

    // Oblivious MUL cost
    static Cost mulCost(bool aShare=true) {
        if (aShare) {
            Cost c = Cost(1,1);
            return c;
        }
        Cost c = b2aCost();
        c.operation++;
        c.synchronization++;
        return c;
    }

    // Oblivious PLUS cost
    static Cost plusCost(bool aShare=true) {
        if (aShare) {
            Cost c = Cost(1,0);
            return c;
        }
        Cost c = b2aCost();
        c.operation++;
        return c;
    }

    // Oblivious RCA cost
    static Cost rcaCost() {
        return Cost(5*LENGTH - 3, LENGTH);
    }

    // Oblivious a2b cost
    static Cost a2bCost() {
        return rcaCost();
    }

    // Oblivious b2a cost
    static Cost b2aCost() {
        return rcaCost();
    }

    // Oblivious single-bit conversion cost
    static Cost bitConvCost(card_t cardinality=1) {
        // NOTE: Single-bit conversion needs roughly the same amount of local operations as logical AND
        double batches = ceil((cardinality) / BATCH_SIZE);
        Cost c = Cost(cardinality * 1, 2 * batches);
        return c;
    }

    // Oblivious equality cost
    static Cost eqCost(card_t cardinality=1) {
        double batches = ceil((cardinality) / BATCH_SIZE);
        Cost c = Cost((2*LENGTH - 1) * cardinality, log2(LENGTH) * batches);
        return c;
    }

    // Oblivious inequality cost
    static Cost ineqCost(card_t cardinality=1) {
        double batches = ceil((cardinality) / BATCH_SIZE);
        Cost c = Cost((4*LENGTH - 3) * cardinality, (log2(LENGTH) + 1) * batches);
        return c;
    }

    // Oblivious cmp-n-swap cost
    static Cost cmpCost(card_t cardinality=1) {
        double batches = ceil(cardinality / BATCH_SIZE);
        Cost c = ineqCost(cardinality);
        c.operation += batches;
        c.synchronization += batches;
        return c;
    }

    // Scan operator costs
    static Cost scanCost(double inCard) {
        return Cost(); // Scan is free
    }

    // Select operator costs
    static Cost selectCost(double inCard, Cost& predCost, unsigned compBits=0) {
        compBits++;
        Cost c = Cost();
        double batches = ceil(inCard / BATCH_SIZE);
        c.operation = inCard * predCost.operation;
        c.synchronization = batches * predCost.synchronization;
        // Additional cost due to composition
        c.compOp = inCard * ceil(log2(compBits));
        c.compSync = batches * ceil(log2(compBits));
        return c;
    }

    // Nested-loop join operator costs
    static Cost joinCost(double lCard, double rCard, Cost& predCost, unsigned compBits=0) {
        compBits++;
        Cost c = Cost();
        double batches = ceil((lCard * rCard) / BATCH_SIZE);
        c.operation = lCard * rCard * predCost.operation;
        c.synchronization = batches * predCost.synchronization;
        // Additional cost due to composition
        c.compOp = lCard * rCard * ceil(log2(compBits));
        c.compSync = batches * ceil(log2(compBits));
        return c;
    }

    // Semi-join operator costs
    static Cost semiJoinCost(double lCard, double rCard, Cost& predCost, unsigned leftCompBits=0,
                             unsigned rightCompBits=0) {
        leftCompBits++;
        rightCompBits++;
        Cost c = Cost();
        double batches = ceil((lCard * rCard) / BATCH_SIZE);
        c.operation = (lCard * rCard * predCost.operation) + lCard * (rCard - 1);
        c.synchronization = batches * (predCost.synchronization + ceil(log2(rCard)));
        // Additional cost due to composition (left)
        c.compOp = lCard * ceil(log2(leftCompBits));
        c.compSync = ceil(lCard / BATCH_SIZE) * ceil(log2(leftCompBits));
        // Additional cost due to composition (right)
        c.compOp += (rCard - 1) * ceil(log2(leftCompBits));
        c.compSync += batches * ceil(log2(rightCompBits));
        return c;
    }

    // Oblivious sort cost
    static Cost sortCost(double inCard) {
        Cost c = Cost(), inqCost = ineqCost();
        double batches = ceil((inCard / 2) / BATCH_SIZE);
        c.operation = 0.25 * inCard * log2(inCard) * (log2(inCard) + 1) * (inqCost.operation + 6);
        c.synchronization = 2 * log2(inCard) * (log2(inCard) + 1) * inqCost.synchronization * batches;
        return c;
    }

    // Group-by operator costs
    static Cost groupByCost(double inCard, Cost& aggCost, unsigned compBits=0) {
        compBits++;
        Cost c = Cost();
        Cost sCost = sortCost(inCard);
        double batches = ceil((inCard * (log2(inCard) - 1) + 1) / BATCH_SIZE);
        c.operation = sCost.operation + inCard * ((log2(inCard) - 1) + 1) * aggCost.operation;
        c.synchronization = sCost.synchronization + log2(inCard) * aggCost.synchronization * batches;
        // Additional cost due to composition
        c.compOp = (inCard * (log2(inCard) - 1) + 1) * ceil(log2(compBits));
        c.compSync = batches * ceil(log2(compBits));
        return c;
    }

    // Odd-even merge costs
    static Cost oddEvenMerge(double inCard, Cost& aggCost, unsigned compBits=0) {
        compBits++;
        Cost c = Cost();
        double batches = ceil((inCard * (log2(inCard) - 1) + 1) / BATCH_SIZE);
        c.operation = inCard * ((log2(inCard) - 1) + 1) * aggCost.operation;
        c.synchronization = log2(inCard) * aggCost.synchronization * batches;
        // Additional cost due to composition
        c.compOp = (inCard * (log2(inCard) - 1) + 1) * ceil(log2(compBits));
        c.compSync = batches * ceil(log2(compBits));
        return c;
    }

    // Distinct operator costs
    static Cost distinctCost(unsigned inCard, unsigned compBits=0) {
        compBits++;
        Cost c = Cost(), sCost = sortCost(inCard), eCost = eqCost();
        double batches = ceil(inCard / BATCH_SIZE);
        c.operation = sCost.operation + (inCard - 1) * eCost.operation;
        c.synchronization = sCost.synchronization + (eCost.synchronization * batches);
        // Masking cost
        c += masking(inCard);
        // Additional cost due to composition
        c.compOp = inCard * ceil(log2(compBits));
        c.compSync = batches * ceil(log2(compBits));
        return c;
    }

    // Distinct adjacent eq costs
    static Cost masking(unsigned inCard, unsigned cols=0) {
        return eqCost(inCard);
    }

    // Distinct adjacent eq costs
    static Cost adjacentEq(unsigned inCard, unsigned compBits=0) {
        compBits++;
        Cost c = Cost(), eCost = eqCost();
        double batches = ceil(inCard / BATCH_SIZE);
        c.operation = (inCard - 1) * eCost.operation;
        c.synchronization = eCost.synchronization * batches;
        // Additional cost due to composition
        c.compOp = inCard * ceil(log2(compBits));
        c.compSync = batches * ceil(log2(compBits));
        return c;
    }

    // Distinct adjacent eq costs
    static Cost adjacentGeq(unsigned inCard, unsigned compBits=0) {
        compBits++;
        Cost c = Cost(), geqCost = ineqCost();
        double batches = ceil(inCard / BATCH_SIZE);
        c.operation = (inCard - 1) * geqCost.operation;
        c.synchronization = geqCost.synchronization * batches;
        // Additional cost due to composition
        c.compOp = inCard * ceil(log2(compBits));
        c.compSync = batches * ceil(log2(compBits));
        return c;
    }
};
#endif  // SECRECY_COSTMODEL_H