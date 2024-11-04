#include "vector.hpp"

#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"

std::pair<size_t, size_t> findCoordsOfFirstLeaf(const Ancestry &ancestry,
                                                int leaf) {
    std::pair<size_t, size_t> coords;
    for (size_t r = 0; r < ancestry.size(); ++r) {
        for (size_t c = 0; c < 3; ++c) {
            if (ancestry[r][c] == leaf) {
                coords = std::make_pair(r, c);
                break;
            }
        }
    }

    return coords;
}

void addLeaf(PhyloVec &v, unsigned int leaf, unsigned int pos) {
    // Append a new leaf branching out from "pos"
    v.push_back(pos);

    Ancestry ancestryAdd = getAncestry(v);

    // Block l. 23 to l. 44 should be simplifiable
    std::pair<size_t, size_t> leafCoords =
        findCoordsOfFirstLeaf(ancestryAdd, v.size());

    size_t leafRow = leafCoords.first;
    size_t leafCol = leafCoords.second;

    // Set the leaf value in the ancestry to - 1
    // (Temporary masking)
    ancestryAdd[leafRow][leafCol] = -1;

    // Increment the other leaves that are >= leaf
    for (size_t r = 0; r < ancestryAdd.size(); ++r) {
        for (size_t c = 0; c < 3; ++c) {
            if (ancestryAdd[r][c] >= leaf) {
                ancestryAdd[r][c] += 1;
            }
        }
    }

    // Re-instate the missing leaf
    ancestryAdd[leafRow][leafCol] = leaf;

    // Re-order the cherries
    orderCherries(ancestryAdd);

    orderCherriesNoParents(ancestryAdd);

    // Build the new vector
    v = buildVector(ancestryAdd);
}

unsigned int removeLeaf(PhyloVec &v, unsigned int leaf) {
    Ancestry ancestry = getAncestry(v);

    std::pair<size_t, size_t> leafCoords =
        findCoordsOfFirstLeaf(ancestry, leaf);

    size_t leafRow = leafCoords.first;
    size_t leafCol = leafCoords.second;

    // Find the parent of the leaf to remove
    unsigned int parent = ancestry[leafRow][2];
    unsigned int sister = ancestry[leafRow][1 - leafCol];

    const size_t numCherries = ancestry.size();
    Ancestry ancestryRm(numCherries - 1);

    for (size_t r = 0; r < numCherries - 1; ++r) {
        if (r < leafRow) {
            ancestryRm[r] = ancestry[r];
        } else {
            ancestryRm[r] = ancestry[r + 1];
        }

        for (size_t c = 0; c < 3; ++c) {
            unsigned int node = ancestryRm[r][c];
            if (node == parent) {
                node = sister;
            }

            // Subtract 1 for leaves > "leaf"
            // (so that the vector is still valid)
            if (node > leaf) {
                node -= 1;

                if (node >= parent) {
                    node -= 1;
                }
            }

            ancestryRm[r][c] = node;
        }
    }

    // We now have a valid ancestry without "leaf"
    // So we build a vector from it

    orderCherries(ancestryRm);

    orderCherriesNoParents(ancestryRm);

    v = buildVector(ancestryRm);

    return sister;
}

std::vector<std::vector<unsigned int>> getAncestryPaths(const PhyloVec &v) {
    Ancestry ancestry = getAncestry(v);

    unsigned int root = 2 * v.size();

    std::vector<unsigned int> parent_vec(root, 0);

    for (size_t i = 0; i < v.size(); ++i) {
        parent_vec[ancestry[i][0]] = ancestry[i][2];
        parent_vec[ancestry[i][1]] = ancestry[i][2];
    }

    std::vector<std::vector<unsigned int>> ancestryPaths;

    for (unsigned int j = 0; j < root; ++j) {
        unsigned int lastNode = j;
        std::vector<unsigned int> path = {lastNode};

        while (lastNode != root) {
            lastNode = parent_vec[path[0]];
            path.push_back(lastNode);
        }

        ancestryPaths.push_back(path);
    }

    return ancestryPaths;
}
int getCommonAncestor(const PhyloVec &v, unsigned int node1,
                      unsigned int node2) {
    std::vector<std::vector<unsigned int>> ancestry_paths = getAncestryPaths(v);

    std::vector<unsigned int> path1 = ancestry_paths[node1];
    std::vector<unsigned int> path2 = ancestry_paths[node2];

    unsigned int i = 0, j = 0;
    while (i < path1.size() && j < path2.size()) {
        if (path1[i] == path2[j]) {
            return path1[i]; // Found MRCA
        } else if (path1[i] < path2[j]) {
            i++; // Move pointer in vec1 forward
        } else {
            j++; // Move pointer in vec2 forward
        }
    }

    return 2 * v.size(); // No common element found --> return root
}