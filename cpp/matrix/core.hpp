#ifndef MATRIX_CORE_HPP
#define MATRIX_CORE_HPP

/**
 * @file core.hpp
 * @brief Core types necessary for Phylo2Mat
 */

#include <array>
#include <vector>

struct PhyloMat {
    PhyloVec v;
    std::vector<std::array<float, 2>> branches;
};

#endif // MATRIX_CORE_HPP