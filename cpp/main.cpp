#include "base/core.hpp"
#include "base/to_newick.hpp"
#include "base/to_vector.hpp"
#include "utils/random.hpp"

#include <iostream>

int main() {
    PhyloVec v = sample(7, false);

    for (auto elm : v) {
        std::cout << elm << " ";
    }

    std::cout << std::endl;

    std::string newick = toNewick(v);

    std::cout << newick << std::endl;

    PhyloVec vec = toVector(newick);

    std::string_view newickNoParents = "((0,5),((1,(3,4)),2));";
    // newickNoParents = "(0,((1,3),(2,4)));";

    // newickNoParents = "(((0,6),2),(((1,5),4),3));";

    Ancestry ancNoParents = reduceNoParents(newickNoParents);

    std::cout << "ancestry" << std::endl;
    for (auto elm : ancNoParents) {
        for (auto sub_elm : elm) {
            std::cout << sub_elm << " ";
        }
        std::cout << std::endl;
    }

    return 0;
}