#include <iostream>

#include "base/core.hpp"
#include "base/to_newick.hpp"
#include "base/to_vector.hpp"
#include "ops/newick.hpp"
#include "ops/vector.hpp"

int main() {
    PhyloVec v = sample(2048, false);

    // std::cout << "v = ";
    // for (auto elm : v) {
    //     std::cout << elm << " ";
    // }

    // std::cout << std::endl;

    std::string newick = toNewick(v);

    removeParentLabels(newick);

    // std::cout << "newick: " << newick << std::endl;

    Ancestry anc = getCherriesNoParents(newick.substr(0, newick.length() - 1));

    orderCherriesNoParents(anc);

    // for (auto row : anc) {
    //     for (auto elm : row) {
    //         std::cout << elm << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << std::endl;

    return 0;
}