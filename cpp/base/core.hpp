#ifndef CORE_HPP
#define CORE_HPP

#include <array>
#include <string_view>
#include <vector>

typedef std::vector<unsigned int> PhyloVec;
typedef std::vector<std::array<int, 3>> Ancestry;
// index = leaf, string = taxon
typedef std::vector<std::string_view> Leaf2Taxon;

#endif // CORE_HPP