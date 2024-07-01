#include <vector>

#include <gtest/gtest.h>

#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"
#include "../utils/random.hpp"
// #include "../utils/validation.hpp"

const int MIN_N_LEAVES = 5;
const int MAX_N_LEAVES = 200;
const int N_REPEATS = 10;

class V2Newick2VTest : public ::testing::TestWithParam<int>
{
protected:
};

INSTANTIATE_TEST_SUITE_P(RandomTests, V2Newick2VTest, ::testing::Range(MIN_N_LEAVES, MAX_N_LEAVES));

TEST_P(V2Newick2VTest, TestV2Newick2V)
{
    int n_leaves = GetParam();
    for (int _ = 0; _ < N_REPEATS; _++)
    {
        PhyloVec v = sample(n_leaves);
        std::string newick = toNewick(v);

        PhyloVec v2 = toVector(newick);

        ASSERT_TRUE(std::equal(v.begin(), v.end(), v2.begin()));
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
