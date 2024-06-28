#ifndef VALIDATION_HPP
#define VALIDATION_HPP

#include "base/core.hpp"

/**
 * @brief check that Phylo2Vec v is correct
 * i.e. that each 0 <= v[i] <= 2i
 */
void check_v(const PhyloVec &v);

#endif // VALIDATION_HPP