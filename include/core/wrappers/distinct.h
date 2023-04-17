#ifndef SECRECY_DISTINCT_H
#define SECRECY_DISTINCT_H

#include "../../common/table.h"
#include "../../common/communication.h"
#include "../relational.h"

// TODO: this is a duplicate of "distinct_batch"
//  - This one is with template.
//  - It has the communication round in the end.
template<typename T, typename T2>
std::shared_ptr<T *> distinct_array(T2 *table, unsigned *key_indices, unsigned num_keys, int batchSize) {

    auto shares = allocateColumnWiseContent<T>(table->numRows, 1);
    auto shares_ptr = shares.get();
    T *distinct = shares_ptr[0];
    int rows = table->numRows;

    assert((batchSize > 0) && (batchSize <= table->numRows - 1));
    BShare **c = table->content;
    BShare mask = 1, max = 0xFFFFFFFFFFFFFFFF;
    // Allocate arrays for local and remote shares of equality results
    BShare *local = (BShare *) malloc(batchSize * num_keys * sizeof(BShare));
    assert(local != NULL);
    BShare *remote = (BShare *) malloc(batchSize * num_keys * sizeof(BShare));
    assert(remote != NULL);


    // First element of the sorted table is always in the set of distinct elements
    distinct[0] = get_rank() / 2;
    int num_bits = sizeof(BShare) * 8;
    int num_levels = log2(num_bits);
    int total_comparisons = table->numRows - 1;
    // For each 'sliding' batch
    int next_start, batch_comp, offset, att_index;
    for (int i = 0; i < total_comparisons; i += batchSize) {
        // The start index of the next batch
        next_start = i + batchSize;
        // The number of comparisons in the current batch
        batch_comp = next_start <= total_comparisons ?
                     batchSize : total_comparisons - i;
        for (int j = 0, k = i; j < batch_comp; j++, k++) {
            offset = j * num_keys;
            for (int idx = 0; idx < num_keys; idx++) {
                att_index = key_indices[idx];
                // x_i ^ y_i ^ 1
                local[offset + idx] = c[k][att_index] ^ c[k + 1][att_index] ^ max;
                remote[offset + idx] = c[k][att_index + 1] ^ c[k + 1][att_index + 1] ^ max;
            }
        }
        // Apply equalities in bulk
        eq_bulk(num_levels, num_bits, local, remote, batch_comp * num_keys);
        exchange_shares_array(local, remote, batch_comp * num_keys);  // 1 round
        // AND all equality bits per pair of adjacent rows to compute the final bit bs
        and_b_all_group(local, remote, batch_comp, num_keys);   // log(num_keys) rounds
        // Set results
        for (int j = 0, k = i + 1; j < batch_comp; j++, k++) {
            distinct[k] = local[j] ^ mask;
        }
    }
    free(local);
    free(remote);

    // communication round
    exchangeArrayWithParties<T>(shares_ptr[0], shares_ptr[1], rows, 1, Distinct_SHARE_TAG, 1);

    return shares;
}

#endif //SECRECY_DISTINCT_H
