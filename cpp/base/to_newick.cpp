#include "to_newick.hpp"

Ancestry getAncestry(const PhyloVec &v) {
    // numLeaves - 1
    const size_t k = v.size();

    // Matrix with 3 columns: child1, child2, parent
    Ancestry ancestry(k);

    std::vector<std::array<unsigned int, 2>> pairs;
    pairs.reserve(k);

    // The first pair of leaves if always (0, 1)
    pairs.push_back({0, 1});

    // The goal here is to add mergers like in the previous iteration
    for (size_t i = 1; i < k; ++i) {
        // The nextLeaf to add is i + 1
        unsigned int nextLeaf = i + 1;

        if (v[i] <= i) {
            /*
            If v[i] <= i, it's like a birth-death process
            We now that the next pair to add now is (v[i], next_leaf)
            (as the branch leading to v[i] gives birth to the next_leaf)

            Why pairs.insert(0)? Let's take an example with [0, 0]
            We initially have (0, 1), but 0 gives birth to 2 afterwards
            So the "shallowest" pair is (0, 2)
            */
            pairs.insert(pairs.begin(), {v[i], nextLeaf});
        } else {
            /*
            If v[i] > i, it's not the branch leading v[i] that gives birth but
            an internal branch

            Thus, it will not be the "shallowest" pair, so we do not insert it
            at position 0, but instead at v[i] - i (i = len(pairs))

            pairs[v[i] - i - 1][0] is a node that we processed
            beforehand which is deeper than the branch v[i]
            */
            pairs.insert(pairs.begin() + v[i] - i,
                         {pairs[v[i] - i - 1][0], nextLeaf});
        }
    }

    // Keep track of the following relationship: child->highest parent
    std::vector<int> parents(2 * k + 1, -1);

    for (size_t i = 0; i < k; ++i) {
        auto &[c1, c2] = pairs[i];

        int parentOfChild1 = parents[c1] != -1 ? parents[c1] : c1;
        int parentofChild2 = parents[c2] != -1 ? parents[c2] : c2;

        // Next parent
        int nextParent = k + i + 1;
        ancestry[i] = {parentOfChild1, parentofChild2, nextParent};

        // Change the parents of the current children
        parents[c1] = nextParent;
        parents[c2] = nextParent;
    }

    return ancestry;
}

static std::string buildNewickRecursiveInner(int p, const Ancestry &ancestry) {
    const int leafMax = ancestry.size();

    auto &[c1, c2, _] = ancestry[p - leafMax - 1];

    std::string left = c1 > leafMax ? buildNewickRecursiveInner(c1, ancestry)
                                    : std::to_string(c1);
    std::string right = c2 > leafMax ? buildNewickRecursiveInner(c2, ancestry)
                                     : std::to_string(c2);

    std::string newick = "(" + left + "," + right + ")" + std::to_string(p);

    return newick;
}

std::string buildNewick(const Ancestry &ancestry) {
    const int root = ancestry.back()[2];

    return buildNewickRecursiveInner(root, ancestry) + ";";
}

std::string toNewick(const PhyloVec &v) {
    Ancestry anc = getAncestry(v);
    return buildNewick(anc);
}