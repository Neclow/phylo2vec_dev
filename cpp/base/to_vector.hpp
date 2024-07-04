#ifndef TO_VECTOR_HPP
#define TO_VECTORCORE_HPP

#include <array>
#include <string>
#include <string_view>
#include <vector>

#include "core.hpp"

/**
 * @brief get all "cherries" from a Newick string with parents
 *
 * @code
 * std::string newick = (((0,4)5,2)7,(1,3)6)8;
 * Ancestry cherries = getCherries(newick);
 * // cherries =
 * // 0 4 5
 * // 5 2 7
 * // 1 3 6
 * // 7 6 8
 * @endcode
 * @param newick Newick string with parents
 * @return Ancestry: vector of triplets {child1, child1, parent}
 */
Ancestry getCherries(std::string_view newick);

Ancestry getCherriesNoParents(std::string_view newick);

/**
 * @brief order all cherries according to their height
 * (i.e., from leaf-level cherries to the root-level pairs)
 * while traversing the cherry-list, we replace the internal node
 * by their smallest descending leaf
 *
 * @code
 * // cherries =
 * // 0 4 5
 * // 5 2 7
 * // 1 3 6
 * // 7 6 8
 * orderCherries(cherries);
 * // output =
 * // 0 4 4
 * // 1 3 3
 * // 0 2 2
 * // 0 1 1
 * @endcode
 * @param ancestry
 */
void orderCherries(Ancestry &ancestry);

void orderCherriesNoParents(Ancestry &ancestry);

PhyloVec buildVector(const Ancestry &cherries);

PhyloVec toVector(std::string_view newick);

PhyloVec toVectorNoParents(std::string_view newick_no_parents);

#endif // TO_VECTOR_HPP