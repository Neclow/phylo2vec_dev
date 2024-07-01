#include "random.hpp"
#include <random>

PhyloVec sample(const int &n_leaves, bool ordered)
{
    PhyloVec v;
    std::random_device rd;
    std::mt19937 gen(rd());

    if (ordered)
    {
        for (std::size_t i = 0; i < n_leaves - 1; ++i)
        {
            std::uniform_int_distribution<> distrib(0, 2 * i);
            v.push_back(distrib(gen));
        }
    }
    else
    {
        for (std::size_t i = 0; i < n_leaves - 1; ++i)
        {
            std::uniform_int_distribution<> distrib(0, i);
            v.push_back(distrib(gen));
        }
    }

    return v;
}