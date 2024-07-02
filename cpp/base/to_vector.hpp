#ifndef TO_VECTOR_HPP
#define TO_VECTORCORE_HPP

#include <array>
#include <string>
#include <string_view>
#include <vector>

#include "core.hpp"

size_t getNumLeavesFromNewick(std::string_view newick);

void doReduce(Ancestry &ancestry, std::string &newick);

Ancestry reduce(std::string &newick);

Ancestry reduceNoParents(std::string newick);

void toCherries(Ancestry &ancestry);

Ancestry orderCherriesNoParents(Ancestry cherries);

PhyloVec buildVector(const Ancestry &cherries);

PhyloVec toVector(std::string &newick);

PhyloVec toVectorNoParents(std::string &newick_no_parents);

#endif // TO_VECTOR_HPP