#include "to_newick.hpp"

Ancestry getAncestry(const PhyloVec &v) {
    // numLeaves - 1
    size_t k = v.size();

    // Matrix with 3 columns: child1, child2, parent
    Ancestry ancestry(k);

    std::vector<std::array<unsigned int, 2>> pairs;
    pairs.reserve(k);

    // The first pair of leaves if always (0, 1)
    pairs.push_back({0, 1});

    // The nextLeaf to add is always 2
    unsigned int nextLeaf = 2;

    // The goal here is to add mergers like in the previous iteration
    for (std::size_t i = 1; i < k; ++i) {
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

        // Increment nextLeaf
        ++nextLeaf;
    }

    // Keep track of the following relationship: child->highest parent
    std::vector<int> parents(k * 2 + 1, -1);
    // Keep track of siblings (i.e., sister nodes)
    // std::vector<int> siblings(v.size() * 2 + 1, -1);

    // next parent
    int nextParent = k + 1;
    int parentOfChild1, parentofChild2;

    // int siblingOfChild1, siblingofChild2;

    for (std::size_t i = 0; i < k; ++i) {
        auto &[c1, c2] = pairs[i];

        parentOfChild1 = parents[c1] != -1 ? parents[c1] : c1;
        parentofChild2 = parents[c2] != -1 ? parents[c2] : c2;

        ancestry[i] = {parentOfChild1, parentofChild2, nextParent};

        // Change the parents of the current children
        parents[c1] = nextParent;
        parents[c2] = nextParent;

        // TODO: this whole block seems irrelevant? All tests pass without it

        // siblingOfChild1 = siblings[c1] != -1 ? siblings[c1] : c1;
        // siblingofChild2 = siblings[c2] != -1 ? siblings[c2] : c2;

        // // Change the parents of the siblings
        // parents[siblingOfChild1] = nextParent;
        // parents[siblingofChild2] = nextParent;

        // // Change the previous parents of the child if there are any
        // parents[parentOfChild1] = nextParent;
        // parents[parentofChild2] = nextParent;

        // // Update siblings
        // siblings[c1] = c2;
        // siblings[c2] = c1;

        ++nextParent;
    }

    return ancestry;
}

std::string buildNewick(const Ancestry &ancestry) {
    // Row with 2 children of root + root node
    auto &[c1, c2, p] = ancestry.back();

    // Convert the ints to strings
    std::string c1Str = std::to_string(c1);
    std::string c2Str = std::to_string(c2);
    std::string pStr = std::to_string(p);

    // Initial Newick: (c1, c2)p;
    std::string newick = "(" + c1Str + "," + c2Str + ")" + pStr + ";";

    // Max leaf number
    int leafMax = ancestry.size();

    // Queue of internal nodes
    std::vector<int> parentQueue;
    // Push nodes in the queue if they are internal
    if (c1 > leafMax) {
        parentQueue.push_back(c1);
    }
    if (c2 > leafMax) {
        parentQueue.push_back(c2);
    }

    // position of nodes in the string
    std::vector<int> nodeIdxs(2 * leafMax);
    nodeIdxs[c1] = 1;
    nodeIdxs[c2] = 2 + c1Str.length();

    // Next parental node to handle
    int nextParent;

    for (int i = 1; i < leafMax; ++i) {
        // Get the next internal node and pop it from the queue
        nextParent = parentQueue.back();
        parentQueue.pop_back();

        // Get the row corresponding to the nextParent
        // TODO: I think p == nextParent
        auto &[c1, c2, p] = ancestry[nextParent - leafMax - 1];

        // Convert the ints to strings
        c1Str = std::to_string(c1);
        c2Str = std::to_string(c2);
        pStr = std::to_string(p);

        // Insert a sub-newick (c1, c2) where the parent is
        // Bef insert: old newick = (0,7)8; c1 = 6, c2 = 2, p = 7
        newick.insert(nodeIdxs[p], "(" + c1Str + "," + c2Str + ")");
        // Aft insert : new newick = (0,(6,2)7)8;

        // Update position of c1, c2 in the string
        nodeIdxs[c1] = nodeIdxs[p] + 1;
        nodeIdxs[c2] = nodeIdxs[c1] + c1Str.length() + 1;

        // Push nodes in the queue if they are internal
        if (c1 > leafMax) {
            parentQueue.push_back(c1);
        }
        if (c2 > leafMax) {
            parentQueue.push_back(c2);
        }
    }

    return newick;
}

std::string toNewick(const PhyloVec &v) { return buildNewick(getAncestry(v)); }