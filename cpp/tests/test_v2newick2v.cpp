#include <gtest/gtest.h>

#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"
#include "../ops/newick.hpp"
#include "../ops/vector.hpp"
#include "config.cpp"

class V2Newick2VTest : public ::testing::TestWithParam<int> {
   protected:
};

INSTANTIATE_TEST_SUITE_P(RandomTests, V2Newick2VTest, ::testing::Range(MIN_N_LEAVES, MAX_N_LEAVES));

TEST_P(V2Newick2VTest, V2Newick2V) {
    int numLeaves = GetParam();
    for (size_t _ = 0; _ < N_REPEATS; ++_) {
        PhyloVec v = sample(numLeaves, false);
        std::string newick = toNewick(v);

        PhyloVec v2 = toVector(newick);

        EXPECT_EQ(v, v2);
    }
}

TEST_P(V2Newick2VTest, Cherries) {
    int numLeaves = GetParam();
    PhyloVec v = sample(numLeaves, false);
    std::string newick = toNewick(v);

    Ancestry anc = getCherries(newick);

    orderCherries(anc);

    removeParentLabels(newick);

    Ancestry ancNoParents = getCherriesNoParents(newick);

    orderCherriesNoParents(ancNoParents);

    EXPECT_EQ(anc, ancNoParents);
}
