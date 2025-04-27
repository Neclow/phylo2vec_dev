#include "to_newick.hpp"

AVLTree makeTree(const PhyloVec &v) {
    const size_t k = v.size();

    AVLTree avl_tree;

    avl_tree.insert(0, {0, 1});

    for (size_t i = 1; i < k; ++i) {
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
            avl_tree.insert(0, {v[i], nextLeaf});
        } else {
            /*
            If v[i] > i, it's not the branch leading v[i] that gives birth
            but an internal branch

            Thus, it will not be the "shallowest" pair, so we do not insert
            it at position 0, but instead at v[i] - i (i = len(pairs))

            pairs[v[i] - i - 1][0] is a node that we processed
            beforehand which is deeper than the branch v[i]
            */
            unsigned int index = v[i] - nextLeaf;
            Pair value = avl_tree.lookup(avl_tree.getRoot(), index);

            avl_tree.insert(index + 1, {value[0], nextLeaf});
        }
    }

    return avl_tree;
}

Pairs getPairs(const PhyloVec &v) {
    AVLTree tree = makeTree(v);
    return tree.getPairs();
}

[[deprecated("getAncestry is no longer used in toNewick, and is left for legacy reasons")]] Ancestry
getAncestry(const PhyloVec &v) {
    const size_t k = v.size();

    Pairs pairs = getPairs(v);

    // Matrix with 3 columns: child1, child2, parent
    Ancestry ancestry(k);

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

std::string buildNewick(const Pairs &pairs) {
    const unsigned int numLeaves = pairs.size() + 1;

    std::vector<std::string> cache;
    cache.reserve(numLeaves);
    for (size_t i = 0; i < numLeaves; ++i) {
        cache.push_back(std::to_string(i));
    }

    for (size_t i = 0; i < pairs.size(); ++i) {
        auto &[c1, c2] = pairs[i];

        std::string next_parent = std::to_string(numLeaves + i);

        cache[c1] = "(" + std::move(cache[c1]) + "," + std::move(cache[c2]) + ")" + next_parent;
    }

    return cache[0] + ";";
}

std::string toNewick(const PhyloVec &v) {
    Pairs pairs = getPairs(v);
    return buildNewick(pairs);
}