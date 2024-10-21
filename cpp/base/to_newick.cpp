#include "to_newick.hpp"

std::vector<Pair> getPairs(const PhyloVec &v) {
    // numLeaves - 1
    const size_t k = v.size();

    std::vector<Pair> pairs;
    pairs.reserve(k);

    for (unsigned int i = k; i-- > 0;) {
        /*
        If v[i] <= i, it's like a birth-death process
        We now that the next pair to add now is (v[i], next_leaf)
        (as the branch leading to v[i] gives birth to the next_leaf)

        Why pairs.insert(0)? Let's take an example with [0, 0]
        We initially have (0, 1), but 0 gives birth to 2 afterwards
        So the "shallowest" pair is (0, 2)
        */
        unsigned int nextLeaf = i + 1;
        if (v[i] <= i) {
            pairs.push_back({v[i], nextLeaf});
        }
    }

    for (size_t j = 1; j < k; ++j) {
        unsigned int nextLeaf = j + 1;
        if (v[j] == 2 * j) {
            // 2*j = extra root ==> pairing = (0, next_leaf)
            pairs.push_back({0, nextLeaf});
        } else if (v[j] > j) {
            /*
            If v[i] > i, it's not the branch leading v[i] that gives birth
            but an internal branch

            Thus, it will not be the "shallowest" pair, so we do not insert
            it at position 0, but instead at v[i] - i (i = len(pairs))

            pairs[v[i] - i - 1][0] is a node that we processed
            beforehand which is deeper than the branch v[i]
            */
            size_t index = pairs.size() + v[j] - 2 * j;
            pairs.insert(pairs.begin() + index,
                         {pairs[index - 1][0], nextLeaf});
        }
    }

    return pairs;
}

std::vector<Pair> getPairs2(const PhyloVec &v) {
    const size_t k = v.size();

    AVLTree avl_tree;

    avl_tree.insert(0, {0, 1});

    for (size_t i = 1; i < k; ++i) {
        unsigned int nextLeaf = i + 1;

        if (v[i] <= i) {
            avl_tree.insert(0, {v[i], nextLeaf});
        } else {
            unsigned int index = v[i] - nextLeaf;
            Pair value = avl_tree.lookup(avl_tree.root, index);

            avl_tree.insert(index + 1, {value[0], nextLeaf});
        }
    }

    return avl_tree.getPairs();
}

Ancestry getAncestry(const PhyloVec &v) {
    const size_t k = v.size();

    std::vector<Pair> pairs = getPairs2(v);

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