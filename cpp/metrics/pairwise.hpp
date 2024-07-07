#include "../base/core.hpp"

// TODO: this should also cover double-valued matrices
typedef std::vector<std::vector<float>> Matrix;

Matrix copheneticDistances(const PhyloVec &v, bool unrooted = false);

Matrix pairwiseDistances(const PhyloVec &v, std::string_view metric,
                         bool unrooted = false);