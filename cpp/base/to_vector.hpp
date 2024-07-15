#ifndef TO_VECTOR_HPP
#define TO_VECTORCORE_HPP

/**
 * @file to_vector.hpp
 * @brief Newick-to-vector conversion functions
 */

#include "core.hpp"

/**
 * @brief Get all "cherries" from a Newick string with parents
 *
 * Example:
 * ```std::string_view newick = (((0,4)5,2)7,(1,3)6)8;```
 * ```Ancestry cherries = getCherries(newick);```
 * cherries =
 * 0 4 5
 * 5 2 7
 * 1 3 6
 * 7 6 8
 * @param newick Newick string with parents
 * @return Ancestry: vector of triplets {child1, child1, parent}
 */
Ancestry getCherries(std::string_view newick);

/**
 * @brief Get all "cherries" from a Newick string without
 * internal node annotations
 * Example:
 * ```std::string_view newick = (((0,4)5,2)7,(1,3)6)8;```
 * ```Ancestry cherries = getCherriesNoParents(newick);```
 * cherries =
 * 0 4 4
 * 0 2 2
 * 1 3 3
 * 0 1 1
 * @param newick Newick string without parent annotations
 * @return Ancestry: vector of cherry triplets {child1, child2, max(child1,
 * child2)}
 */
Ancestry getCherriesNoParents(std::string_view newick);

/**
 * @brief Order all cherries according to their height
 * (i.e., from leaf-level cherries to the root-level pairs).
 * While traversing the cherry-list, we replace the parent
 * by their largest descending leaf
 *
 * Example:
 * cherries =
 * 0 4 5
 * 5 2 7
 * 1 3 6
 * 7 6 8
 * ```orderCherries(cherries);```
 * output =
 * 0 4 4
 * 1 3 3
 * 0 2 2
 * 0 1 1
 * @param ancestry vector of cherry triplets {child1, child2, max(child1,
 * child2)}
 */
void orderCherries(Ancestry &ancestry);

/**
 * @brief Order all cherries according to their height
 * without internal node annotations
 * (i.e., from leaf-level cherries to the root-level pairs):
 * Example:
 * cherries =
 * 0 4 4
 * 0 2 2
 * 1 3 3
 * 0 1 1
 * ```orderCherriesNoParents(cherries);```
 * output =
 * 0 4 4
 * 1 3 3
 * 0 2 2
 * 0 1 1
 * @param ancestry vector of cherry triplets {child1, child2, max(child1,
 * child2)}
 */
void orderCherriesNoParents(Ancestry &ancestry);

/**
 * @brief construct a Phylo2Vec vector from cherry triplets
 *
 * @param cherries vector of cherry triplets {child1, child2, max(child1,
 * child2)}
 * @return PhyloVec: v[i] = j <=> leaf j descends from branch i
 */
PhyloVec buildVector(const Ancestry &cherries);

/**
 * @brief Convert a newick (with parent annotations) to a Phylo2Vec vector
 * Wrapper of getCherries + orderCherries + buildVector
 * @param newick Newick string with parent labels
 * @return PhyloVec: v[i] = j <=> leaf j descends from branch i
 */
PhyloVec toVector(std::string_view newick);

/**
 * @brief Convert a newick (without parent annotations) to a Phylo2Vec vector
 *Wrapper of getCherries + orderCherriesNoParents + buildVector
 * @param newick_no_parents Newick string without parent labels
 * @return PhyloVec: v[i] = j <=> leaf j descends from branch i
 */
PhyloVec toVectorNoParents(std::string_view newick_no_parents);

#endif // TO_VECTOR_HPP