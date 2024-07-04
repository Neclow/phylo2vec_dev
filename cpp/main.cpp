#include "base/core.hpp"
#include "base/to_newick.hpp"
#include "base/to_vector.hpp"
#include "utils/newick.hpp"
#include "utils/random.hpp"

#include <iostream>

int main() {
    PhyloVec v = sample(2000, false);

    std::string newick = toNewick(v);

    newick = "(((((((tip_0,tip_1),(tip_2,tip_3)),(tip_4,tip_5)),tip_6),tip_7),("
             "((tip_8,tip_9),tip_10),tip_11)),tip_12);";

    std::cout << newick << std::endl;

    Converter conv = toIntNewick(newick);

    std::cout << conv.intNewick << std::endl;

    for (auto elm : conv.mapping) {
        std::cout << elm << " ";
    }
    std::cout << std::endl;

    std::string newick2 = toStringNewick(conv);

    std::cout << newick2 << std::endl;

    return 0;
}