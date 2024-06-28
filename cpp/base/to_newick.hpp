#ifndef TO_NEWICK_HPP
#define TO_NEWICK_HPP

#include <array>
#include <string>
#include <vector>

/**
 * @brief Get ancestry for each node given a v-representation.
 *
 * @param v Phylo2Vec vector
 * @return std::vector<std::array<int, 3>>
 * Ancestry matrix
 * 1st column: child 1
 * 2nd column: child 2
 * 3rd column: parent node
 */
std::vector<std::array<int, 3>> getAncestry(const std::vector<int> &v);

/**
 * @brief
 * Build a Newick string from an "ancestry" array
 * The matrix is processed such that we iteratively write a Newick string
 * to describe the tree.
 * @param ancestry "Ancestry" array of size (n_leaves - 1, 3)
 * @return std::string Newick string
 */
std::string buildNewick(const std::vector<std::array<int, 3>> &ancestry);

/**
 * @brief Wrapper of getAncestry and buildNewick
 *
 * @param v Phylo2Vec vector
 */
std::string toNewick(const std::vector<int> &v);

#endif // TO_NEWICK_HPP