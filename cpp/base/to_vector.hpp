#ifndef TO_VECTOR_HPP
#define TO_VECTOR_HPP

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
 * ===== start =====
 * Iteration 1:
 * 0 4 5 --> 0 4 4
 * child_min[5] = 0
 * Iteration 2:
 * 1 3 6 --> 1 3 3
 * child_min[6] = 1
 * Iteration 3:
 * 5 2 7 --> 0 2 2 because child_min[5] = 0
 * child_min[7] = 0
 * Iteration 4:
 * 7 6 8 --> 0 1 1 because child_min[7] = 0 and child_min[6] = 1
 * child_min[8] = 0
 * ===== end =====
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
 * ===== start =====
 * Iteration 1:
 * 0-4 is the cherry with the highest leaf (4)
 * and both 0 and 4 did not appear previously in the ancestry
 * (processed from top to bottom)
 * output[0] = 0 4 4
 * Iteration 2:
 * Row 0 was processed to it is skipped
 * 1-3 is the cherry with the next highest leaf (3)
 * and both 1 and 3 did not appear previously in the ancestry
 * output[1] = 1 3 3
 * Iteration 3:
 * Rows 0, 2 were processed to they are skipped
 * 0-2 is the cherry with the next highest leaf (2)
 * and both 0 and 2 did not appear previously in the ancestry
 * (as row 0 was skipped)
 * output[2] = 0 2 2
 * Iteration 4:
 * The last row is 0 1 1
 * ===== end =====
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
PhyloVec buildVector(Ancestry cherries);

/**
 * @brief Convert a newick (with parent annotations) to a Phylo2Vec vector
 * Wrapper of getCherries + orderCherries + buildVector
 * @param newick Newick string with parent labels
 * @return PhyloVec: v[i] = j <=> leaf j descends from branch i
 */
PhyloVec toVector(std::string_view newick);

/**
 * @brief Convert a newick (without parent annotations) to a Phylo2Vec vector
 *Wrapper of getCherriesNoParents + orderCherriesNoParents + buildVector
 * @param newick_no_parents Newick string without parent labels
 * @return PhyloVec: v[i] = j <=> leaf j descends from branch i
 */
PhyloVec toVectorNoParents(std::string_view newick_no_parents);

#endif  // TO_VECTOR_HPP