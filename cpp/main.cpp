#include "base/core.hpp"
#include "base/to_newick.hpp"
#include "base/to_vector.hpp"
#include "utils/newick.hpp"
#include "utils/random.hpp"

#include <iostream>

int main() {
    PhyloVec v = sample(2000, false);

    std::string newick = toNewick(v);

    return 0;
}