#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"
#include "../utils/newick.hpp"
#include "../utils/random.hpp"
#include "config.cpp"

#include <gtest/gtest.h>

#include <vector>

class V2Newick2VTest : public ::testing::TestWithParam<int> {
  protected:
};

INSTANTIATE_TEST_SUITE_P(RandomTests, V2Newick2VTest,
                         ::testing::Range(MIN_N_LEAVES, MAX_N_LEAVES));

TEST_P(V2Newick2VTest, V2Newick2V) {
    int numLeaves = GetParam();
    for (size_t _ = 0; _ < N_REPEATS; _++) {
        PhyloVec v = sample(numLeaves, false);
        std::string newick = toNewick(v);

        PhyloVec v2 = toVector(newick);

        ASSERT_TRUE(std::equal(v.begin(), v.end(), v2.begin()));
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

    ASSERT_TRUE(std::equal(anc.begin(), anc.end(), ancNoParents.begin()));
}
