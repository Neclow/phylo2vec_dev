#include "random.hpp"
#include <random>

PhyloVec sample(const size_t &numLeaves, bool ordered) {
    PhyloVec v;
    std::random_device rd;
    std::mt19937 gen(rd());

    if (ordered) {
        for (std::size_t i = 0; i < numLeaves - 1; ++i) {
            std::uniform_int_distribution<> distrib(0, 2 * i);
            v.push_back(distrib(gen));
        }
    } else {
        for (std::size_t i = 0; i < numLeaves - 1; ++i) {
            std::uniform_int_distribution<> distrib(0, i);
            v.push_back(distrib(gen));
        }
    }

    return v;
}