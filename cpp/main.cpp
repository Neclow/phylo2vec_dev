#include "base/core.hpp"
#include "base/to_newick.hpp"
#include "base/to_vector.hpp"
#include "utils/newick.hpp"
#include "utils/random.hpp"

#include <iostream>

int main() {
    // PhyloVec v = sample(5, false);

    // std::string_view newick = toNewick(v);

    std::string_view newick = "(((0,4)5,2)7,(1,3)6)8;";

    PhyloVec v_out = toVector(newick);

    std::cout << "v = ";
    for (auto elm : v_out) {
        std::cout << elm << " ";
    }

    std::cout << std::endl;

    std::string_view newick_no_parents = "(((0,4),2),(1,3));";

    // std::string_view tmp =
    //     "(((((((0,1),(2,3)),(4,5)),6),7),((8,9),10),11)),12);";

    PhyloVec v_out_no_parents = toVectorNoParents(newick_no_parents);

    // // newick = "(((0,4)5,2)7,(1,3)6)8;";

    // newick =
    // "(((((((tip_0,tip_1),(tip_2,tip_3)),(tip_4,tip_5)),tip_6),tip_7),("
    //          "((tip_8,tip_9),tip_10),tip_11)),tip_12);";

    // std::cout << newick << std::endl;

    // Converter conv = toIntNewick(newick);

    // std::cout << conv.intNewick << std::endl;

    // for (auto elm : conv.mapping) {
    //     std::cout << elm << " ";
    // }
    // std::cout << std::endl;

    // std::string newick2 = toStringNewick(conv);

    // std::cout << newick2 << std::endl;

    return 0;
}