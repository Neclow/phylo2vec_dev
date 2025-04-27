#include <gtest/gtest.h>

#include <random>

#include "../ops/vector.hpp"
#include "config.cpp"

class UtilsTest : public ::testing::TestWithParam<int> {
   protected:
};

INSTANTIATE_TEST_SUITE_P(RandomTests, UtilsTest, ::testing::Range(MIN_N_LEAVES, MAX_N_LEAVES));

TEST_P(UtilsTest, SampleTest) {
    int numLeaves = GetParam();
    for (size_t _ = 0; _ < N_REPEATS; ++_) {
        PhyloVec v = sample(numLeaves, false);

        EXPECT_NO_THROW(check_v(v));
    }
}

TEST_P(UtilsTest, RemoveAndAddTest) {
    int numLeaves = GetParam();
    std::random_device rd;
    std::mt19937 gen(rd());
    for (size_t _ = 0; _ < N_REPEATS; ++_) {
        PhyloVec v = sample(numLeaves, false);

        PhyloVec vOld = v;

        std::uniform_int_distribution<> distrib(0, numLeaves - 1);

        unsigned int leaf = distrib(gen);

        unsigned int sister = removeLeaf(v, leaf);

        // To deal with increment when removing a leaf
        if (sister >= leaf) {
            sister -= 1;
        }

        addLeaf(v, leaf, sister);

        EXPECT_EQ(v, vOld);
    }
}