#ifndef CORE_HPP
#define CORE_HPP

/**
 * @file core.hpp
 * @brief Core types necessary for Phylo2Vec
 */

#include <array>
#include <string>
#include <string_view>
#include <vector>

/**
 * @brief Phylo2Vec vector v
 * Let n the number of leaves:
 * v is s.t. 0 <= v[i] <= 2*i for i in [0, n - 1]
 */
typedef std::vector<unsigned int> PhyloVec;
/**
 * @brief Ancestry (cherry list)
 * List of triplets that make a tree.
 * Let c1, c2, p = child1, child2, parent
 * Ancestry is either of the form
 * {c1, c2 p}
 * or
 * {c1, c2, max(c1, c2)}
 */
typedef std::vector<std::array<int, 3>> Ancestry;
/**
 * @brief mapping of integer leaves to a distinct taxon
 * vector index = leaf, string = taxon
 *
 */
typedef std::vector<std::string_view> Leaf2Taxon;

struct Converter {
    Leaf2Taxon mapping;
    std::string intNewick;
};

#endif // CORE_HPP