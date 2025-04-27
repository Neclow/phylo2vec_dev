#ifndef FENWICK_HPP
#define FENWICK_HPP

#include <vector>

class FenwickTree {
   public:
    FenwickTree(unsigned int n);

    unsigned int prefix_sum(unsigned int i);
    void update(unsigned int i, unsigned int delta);

   private:
    unsigned int n_leaves;
    std::vector<unsigned int> data;
};

#endif  // FENWICK_HPP