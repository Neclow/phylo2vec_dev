#ifndef TO_NEWICK_HPP
#define TO_NEWICK_HPP

/**
 * @file to_vector.hpp
 * @brief Vector-to-Newick conversion functions
 */

#include "../utils/avl.hpp"
#include "core.hpp"

AVLTree makeTree(const PhyloVec &v);

Pairs getPairs(const PhyloVec &v);

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
Ancestry getAncestry(const PhyloVec &v);

/**
 * @brief
 * Build a Newick string from an "ancestry" array
 * The matrix is processed such that we iteratively write a Newick string
 * to describe the tree.
 * @param ancestry "Ancestry" array of size (n_leaves - 1, 3)
 * @return std::string Newick string
 */
std::string buildNewick(const Ancestry &ancestry);

/**
 * @brief Convert a Phylo2Vec vector to a Newick string.
 * Wraps getAncestry and buildNewick
 *
 *
 * @param v Phylo2Vec vector
 */
std::string toNewick(const PhyloVec &v);

#endif  // TO_NEWICK_HPP