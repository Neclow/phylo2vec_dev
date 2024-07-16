#ifndef MATRIX_TO_NEWICK_HPP
#define MATRIX_TO_NEWICK_HPP

#include "../base/core.hpp"
#include "core.hpp"

std::string buildNewickWithBranches(const Ancestry &ancestry,
                                    std::vector<std::array<float, 2>> branches);

std::string toNewick(const PhyloMat &m);

#endif // MATRIX_TO_NEWICK_HPP