#include "to_newick.hpp"

#include "../base/to_newick.hpp"

std::string
buildNewickWithBranches(const Ancestry &ancestry,
                        std::vector<std::array<float, 2>> branches) {
    auto &[c1, c2, p] = ancestry.back();
    auto &[b1, b2] = branches.back();

    std::string c1Str = std::to_string(c1) + ":" + std::to_string(b1);
    std::string c2Str = std::to_string(c2) + ":" + std::to_string(b2);
    std::string pStr = std::to_string(p);

    // Initial Newick: (c1:b1, c2:b2)p;
    std::string newick = "(" + c1Str + "," + c2Str + ")" + pStr + ";";

    // Max leaf number
    const int leafMax = ancestry.size();

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

    for (int i = 1; i < leafMax; ++i) {
        // Get the next internal node and pop it from the queue
        int nextParent = parentQueue.back();
        parentQueue.pop_back();

        // Get the row corresponding to the nextParent
        // Note that p == nextParent
        auto &[c1, c2, p] = ancestry[nextParent - leafMax - 1];

        auto &[b1, b2] = branches[nextParent - leafMax - 1];

        // Convert the ints to strings
        std::string c1Str = std::to_string(c1) + ":" + std::to_string(b1);
        std::string c2Str = std::to_string(c2) + ":" + std::to_string(b2);

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

std::string toNewick(const PhyloMat &m) {
    Ancestry ancestry = getAncestry(m.v);

    std::string newick = buildNewickWithBranches(ancestry, m.branches);

    return newick;
}