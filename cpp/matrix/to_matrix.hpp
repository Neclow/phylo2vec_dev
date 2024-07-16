#ifndef TO_MATRIX_HPP
#define TO_MATRIX_HPP

#include "../base/core.hpp"
#include "core.hpp"

std::pair<Ancestry, std::vector<std::array<float, 2>>>
getCherriesAndBranches(std::string_view newick);

std::pair<Ancestry, std::vector<std::array<float, 2>>>
getCherriesAndBranchesNoParents(std::string_view newick);

PhyloMat toMatrix(std::string_view newick);

PhyloMat toMatrixNoParents(std::string_view newick);

#endif // TO_MATRIX_HPP