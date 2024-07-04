#include "base/core.hpp"
#include "base/to_newick.hpp"
#include "base/to_vector.hpp"
#include "utils/newick.hpp"
#include "utils/random.hpp"

#include <iostream>

int main() {
    PhyloVec v = sample(2000, false);

    std::string newick = toNewick(v);

    // Ancestry anc = reduce(newick);

    // toCherries(anc);

    // for (auto row : anc) {
    //     for (auto elm : row) {
    //         std::cout << elm << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << std::endl;

    removeParentLabels(newick);

    // std::cout << newick << std::endl;

    Ancestry ancNoParents = reduceNoParents(newick);

    // toCherriesNoParents(ancNoParents);

    // for (auto row : ancNoParents) {
    //     for (auto elm : row) {
    //         std::cout << elm << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // std::cout << std::endl;

    // std::cout << "to cherries" << std::endl;

    toCherriesNoParents(ancNoParents);

    // for (auto row : ancNoParents) {
    //     for (auto elm : row) {
    //         std::cout << elm << " ";
    //     }
    //     std::cout << std::endl;
    // }

    // for (auto elm : v) {
    //     std::cout << elm << " ";
    // }

    // std::cout << std::endl;

    // std::string newick = toNewick(v);

    // std::cout << newick << std::endl;

    // PhyloVec vec = toVector(newick);

    // std::string_view newickNoParents = "((0,5),((1,(3,4)),2));";
    // // newickNoParents = "(0,((1,3),(2,4)));";

    // // newickNoParents = "(((0,6),2),(((1,5),4),3));";

    // Ancestry ancNoParents = reduceNoParents(newickNoParents);

    // std::cout << "ancestry" << std::endl;
    // for (auto elm : ancNoParents) {
    //     for (auto sub_elm : elm) {
    //         std::cout << sub_elm << " ";
    //     }
    //     std::cout << std::endl;
    // }

    return 0;
}