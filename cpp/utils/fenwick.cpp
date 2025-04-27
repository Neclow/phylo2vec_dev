#include "fenwick.hpp"

FenwickTree::FenwickTree(unsigned int n) : n_leaves(n), data(n + 1, 0) {}

unsigned int FenwickTree::prefix_sum(unsigned int i) {
    unsigned int sum = 0;
    while (i > 0) {
        sum += data[i];
        i -= i & -i;
    }
    return sum;
}

void FenwickTree::update(unsigned int i, unsigned int delta) {
    while (i <= n_leaves) {
        data[i] += delta;
        i += i & -i;
    }
}