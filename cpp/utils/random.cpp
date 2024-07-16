#include "random.hpp"
#include <random>

PhyloVec sample(const size_t &numLeaves, bool ordered) {
    PhyloVec v(numLeaves - 1);
    // std::mt19937 gen(std::random_device{}());
    std::minstd_rand gen(std::random_device{}());

    if (ordered) {
        for (size_t i = 0; i < numLeaves - 1; ++i) {
            std::uniform_int_distribution<> distrib(0, 2 * i);
            v[i] = distrib(gen);
        }
    } else {
        for (size_t i = 0; i < numLeaves - 1; ++i) {
            std::uniform_int_distribution<> distrib(0, i);
            v[i] = distrib(gen);
        }
    }

    return v;
}