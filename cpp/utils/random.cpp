#include "random.hpp"
#include <random>

PhyloVec sample(const size_t &numLeaves, bool ordered) {
    PhyloVec v(numLeaves - 1);
    // static std::mt19937 gen(std::random_device{}());
    static std::minstd_rand gen(std::random_device{}());

    if (ordered) {
        for (size_t i = 1; i < numLeaves - 1; ++i) {
            v[i] = gen() % (i + 1);
        }
    } else {
        for (size_t i = 1; i < numLeaves - 1; ++i) {
            v[i] = gen() % (2 * i + 1);
        }
    }

    return v;
}