#ifndef RANDOM_HPP
#define RANDOM_HPP

#include "../base/core.hpp"

/**
 * @brief Sample a random Phylo2Vec v for n leaves
 *
 * @param n_leaves number of leaves
 */
PhyloVec sample(const int &n_leaves, bool ordered = false);

#endif // RANDOM_HPP