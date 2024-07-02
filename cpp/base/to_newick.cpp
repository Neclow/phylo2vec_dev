#include "to_newick.hpp"

Ancestry getAncestry(const PhyloVec &v) {
    // Matrix with 3 columns: child1, child2, parent
    Ancestry ancestry(v.size());

    // This is the first pair, we start with (0, 1)
    // std::vector<std::pair<int, int>> pairs; seems to be slower
    std::vector<std::array<unsigned int, 2>> pairs;
    pairs.reserve(v.size());

    pairs.push_back({0, 1});

    unsigned int nextLeaf;

    // The goal here is to add mergers like in the previous iteration
    for (std::size_t i = 1; i < v.size(); ++i) {
        nextLeaf = i + 1;

        if (v[i] <= i) {
            /*
            If v[i] <= i, it's an easy BD
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
            an internal branch Remark: it will not be the "shallowest" pair, so
            we do not insert it at position 0 len(pairs) = number of pairings we
            did so far So what v[i] - len(pairs) gives us is the depth of the
            next pairing And pairs[v[i] - len(pairs) - 1][0] is a node that we
            processed beforehand which is deeper than the branch v[i]
            */
            pairs.insert(pairs.begin() + v[i] - pairs.size(),
                         {pairs[v[i] - pairs.size() - 1][0], nextLeaf});
        }
    }

    // Keep track of the following relationship: child->highest parent
    std::vector<int> parents(v.size() * 2 + 1, -1);
    // Keep track of siblings (i.e., sister nodes)
    std::vector<int> siblings(v.size() * 2 + 1, -1);

    int nextParent = v.size() + 1;
    int parentOfChild1, parentofChild2;
    int siblingOfChild1, siblingofChild2;

    for (std::size_t i = 0; i < pairs.size(); ++i) {
        auto &[c1, c2] = pairs[i];

        parentOfChild1 = parents[c1] != -1 ? parents[c1] : c1;
        parentofChild2 = parents[c2] != -1 ? parents[c2] : c2;

        ancestry[i] = {parentOfChild1, parentofChild2, nextParent};

        // Change the parents of the current children
        parents[c1] = nextParent;
        parents[c2] = nextParent;

        siblingOfChild1 = siblings[c1] != -1 ? siblings[c1] : c1;
        siblingofChild2 = siblings[c2] != -1 ? siblings[c2] : c2;

        // Change the parents of the siblings
        parents[siblingOfChild1] = nextParent;
        parents[siblingofChild2] = nextParent;

        // Change the previous parents of the child if there are any
        parents[parentOfChild1] = nextParent;
        parents[parentofChild2] = nextParent;

        // Update siblings
        siblings[c1] = c2;
        siblings[c2] = c1;

        ++nextParent;
    }

    return ancestry;
}

std::string buildNewick(const Ancestry &ancestry) {
    // Row with 2 children of root + root node
    auto &[c1, c2, p] = ancestry.back();

    std::string c1Str = std::to_string(c1);
    std::string c2Str = std::to_string(c2);
    std::string pStr = std::to_string(p);

    std::string newick = "(" + c1Str + "," + c2Str + ")" + pStr + ";";

    std::vector<int> queue;

    int nMax = ancestry.size();

    if (c1 > nMax) {
        queue.push_back(c1);
    }
    if (c2 > nMax) {
        queue.push_back(c2);
    }

    std::vector<int> nodeIdxs(2 * nMax);
    nodeIdxs[c1] = 1;
    nodeIdxs[c2] = 2 + c1Str.length();

    int nextParent;

    for (int i = 1; i < nMax; ++i) {
        nextParent = queue.back();

        queue.pop_back();

        auto &[c1, c2, p] = ancestry[nextParent - nMax - 1];

        c1Str = std::to_string(c1);
        c2Str = std::to_string(c2);
        pStr = std::to_string(p);

        std::string_view sub_newick = "(" + c1Str + "," + c2Str + ")" + pStr;

        newick.insert(nodeIdxs[p], sub_newick);
        newick.erase(nodeIdxs[p] + sub_newick.length(), pStr.length());

        nodeIdxs[c1] = nodeIdxs[p] + 1;
        nodeIdxs[c2] = nodeIdxs[c1] + 1 + c1Str.length();

        if (c1 > nMax) {
            queue.push_back(c1);
        }
        if (c2 > nMax) {
            queue.push_back(c2);
        }
    }

    return newick;
}

std::string toNewick(const PhyloVec &v) { return buildNewick(getAncestry(v)); }