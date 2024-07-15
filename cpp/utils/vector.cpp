#include "vector.hpp"

#include "../base/to_newick.hpp"
#include "../base/to_vector.hpp"

std::pair<size_t, size_t> findCoordsOfFirstLeaf(const Ancestry &ancestry,
                                                const int &leaf) {
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

void addLeaf(PhyloVec &v, const unsigned int &leaf, const unsigned int &pos) {
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

unsigned int removeLeaf(PhyloVec &v, const unsigned int &leaf) {
    Ancestry ancestry = getAncestry(v);
    std::pair<size_t, size_t> leafCoords =
        findCoordsOfFirstLeaf(ancestry, leaf);

    size_t leafRow = leafCoords.first;
    size_t leafCol = leafCoords.second;

    // Find the parent of the leaf to remove
    unsigned int parent = ancestry[leafRow][2];
    unsigned int sister = ancestry[leafRow][1 - leafCol];

    Ancestry ancestryRm(ancestry.size() - 1);

    for (size_t r = 0; r < ancestry.size() - 1; ++r) {
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

            // Subtract 1 from the leafs larger > "leaf"
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

    // We now have a correct ancestry without "leaf"
    // So we build a vector from it

    orderCherries(ancestryRm);

    orderCherriesNoParents(ancestryRm);

    v = buildVector(ancestryRm);

    return sister;
}