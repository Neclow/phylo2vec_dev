#include "to_vector.hpp"

#include <algorithm>
#include <charconv>
#include <stdexcept>

int stoi_substr(std::string_view s, size_t start, size_t *end) {
    int value;
    auto [ptr, ec] =
        std::from_chars(s.data() + start, s.data() + s.size(), value);
    if (ec == std::errc()) {
        *end = ptr - s.data();
        return value;
    } else {
        throw std::logic_error("Bad");
    }
}

Ancestry getCherries(std::string_view newick) {
    Ancestry ancestry;

    // Stack of nodes
    std::vector<int> stack;

    for (size_t i = 0; i < newick.length(); ++i) {
        char c = newick[i];
        if (c == ')') {
            ++i;

            // Pop the children nodes from the stack
            int c2 = stack.back();
            stack.pop_back();

            int c1 = stack.back();
            stack.pop_back();

            // Get the parent node after )
            size_t end;
            int p = stoi_substr(newick, i, &end);
            i = end - 1;

            // Add the triplet (c1, c2, p)
            ancestry.push_back({c1, c2, p});

            // Push the parent node to the stack
            // c1 = p;
            stack.push_back(p);
        } else if (c >= '0' && c <= '9') {
            // Get the next node and push it to the stack
            size_t end;
            int node = stoi_substr(newick, i, &end);
            stack.push_back(node);
            i = end - 1;
        }
    }

    return ancestry;
}

Ancestry getCherriesNoParents(std::string_view newick) {
    Ancestry ancestry;
    ancestry.reserve(newick.length());

    std::vector<int> stack;

    for (size_t i = 0; i < newick.length(); ++i) {
        char c = newick[i];

        if (c == ')') {
            // Pop the children nodes from the stack
            int c2 = stack.back();
            stack.pop_back();

            int c1 = stack.back();
            stack.pop_back();

            // No parent annotation --> store the max leaf
            int cMax = std::max(c1, c2);
            ancestry.push_back({c1, c2, cMax});

            // Push the min leaf to the stack
            int cMin = std::min(c1, c2);
            // c1 = cMin;
            stack.push_back(cMin);
        } else if (c >= '0' && c <= '9') {
            // Get the next leaf and push it to the stack
            size_t end;
            int leaf = stoi_substr(newick, i, &end);
            stack.push_back(leaf);
            i = end - 1;
        }
    }

    return ancestry;
}

void orderCherries(Ancestry &ancestry) {
    const size_t numCherries = ancestry.size();
    const size_t numNodes = 2 * numCherries + 1;

    // To keep track + update of the smallest Descendant of each triplet
    // This allows us to reconstruct the triplets as they should appear
    // using the Phylo2Vec construction
    // Note: the first numLeaves indices are not used
    std::vector<int> minDesc(numNodes, -1);

    // Sort the ancestry by their parent node (ascending order)
    std::qsort(ancestry.data(), numCherries, sizeof(std::array<int, 3>),
               [](const void *a, const void *b) {
                   const auto &arr1 =
                       *static_cast<const std::array<int, 3> *>(a);
                   const auto &arr2 =
                       *static_cast<const std::array<int, 3> *>(b);
                   return arr1[2] - arr2[2];
               });

    for (size_t i = 0; i < numCherries; ++i) {
        auto &[c1, c2, p] = ancestry[i];

        // Get the minimum descendant of c1 and c2 (if they exist)
        // minDesc[child_x] doesn't exist, minDesc_x --> child_x
        int minDesc1 = minDesc[c1] != -1 ? minDesc[c1] : c1;
        int minDesc2 = minDesc[c2] != -1 ? minDesc[c2] : c2;

        // Collect the minimum descendant and allocate it to minDesc[parent]
        int descMin = std::min(minDesc1, minDesc2);
        minDesc[p] = descMin;

        // Instead of the parent, we collect the max node
        int descMax = std::max(minDesc1, minDesc2);
        ancestry[i] = {minDesc1, minDesc2, descMax};
    }
}

void orderCherriesNoParents(Ancestry &ancestry) {
    // numCherries = numLeaves - 1
    const int numCherries = ancestry.size();

    // Copy of current ancestry
    const Ancestry oldAncestry = ancestry;

    // Index ordering for the new ancestry
    std::vector<std::uint8_t> idxs(numCherries, 0);

    for (int i = 0; i < numCherries; ++i) {
        // Find the next index to process:
        // The goal is to find the row with the highest leaf
        // where both leaves were previously un-visited
        // why? If a leaf in a cherry already appeared in the ancestry,
        // it means that leaf was already involved in a shallower cherry
        int idx;

        // Initially, all leaves have not been processed
        std::vector<std::uint8_t> unvisited(numCherries + 1, 1);

        // Temporary max leaf
        int maxLeaf = -1;

        for (int j = 0; j < numCherries; ++j) {
            // Row j was processed --> continue
            if (idxs[j]) {
                continue;
            }

            auto &[c1, c2, cMax] = oldAncestry[j];

            // Are c1 and c2 both un-visited?
            if (unvisited[c1] && unvisited[c2]) {
                // If yes, is max(c1, c2) > maxLeaf?
                if (cMax > maxLeaf) {
                    // If yes again, update maxLeaf
                    maxLeaf = cMax;
                    // index to swap <-- j
                    idx = j;
                }
            }

            // c1 and c2 have been processed
            unvisited[c1] = 0;
            unvisited[c2] = 0;
        }
        // Swap the rows for the new ancestry
        // row idx becomes row i
        ancestry[i] = oldAncestry[idx];

        // Row idx has been processed
        idxs[idx] = 1;
    }
}

PhyloVec buildVector(Ancestry cherries) {
    const size_t numCherries = cherries.size();

    PhyloVec v(numCherries, 0);

    // Note: v[0] is always 0
    // but starting with i = 1 makes some tests fail (weird)

    for (size_t i = 0; i < numCherries; ++i) {
        auto &[c1, c2, cMax] = cherries[i];

        int idx = 0;

        for (int j = i - 1; j >= 0; --j) {
            if (cherries[j][2] <= cMax) {
                ++idx;
            } else {
                break;
            }
        }

        v[cMax - 1] = idx == 0 ? std::min(c1, c2) : cMax - 1 + idx;
    }

    return v;
}

PhyloVec toVector(std::string_view newick) {
    Ancestry ancestry = getCherries(newick.substr(0, newick.length() - 1));

    orderCherries(ancestry);

    PhyloVec v = buildVector(ancestry);

    return v;
}

PhyloVec toVectorNoParents(std::string_view newick) {
    Ancestry ancestry =
        getCherriesNoParents(newick.substr(0, newick.length() - 1));

    orderCherriesNoParents(ancestry);

    PhyloVec v = buildVector(ancestry);

    return v;
}