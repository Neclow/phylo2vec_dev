#ifndef RANDOM_HPP
#define RANDOM_HPP

#include <vector>

/**
 * @brief Sample a random Phylo2Vec v for n leaves
 *
 * @param n_leaves number of leaves
 */
std::vector<int> sample(const int &n_leaves, bool ordered);

#endif // RANDOM_HPP