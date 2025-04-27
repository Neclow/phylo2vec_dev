#include "to_vector.hpp"

#include <algorithm>
#include <charconv>
#include <numeric>
#include <stdexcept>
#include <unordered_map>

#include "../utils/fenwick.hpp"

int stoi_substr(std::string_view s, size_t start, size_t *end) {
    int value;
    auto [ptr, ec] = std::from_chars(s.data() + start, s.data() + s.size(), value);
    if (ec == std::errc()) {
        *end = ptr - s.data();
        return value;
    } else {
        throw std::logic_error("Bad");
    }
}

Ancestry getCherries(std::string_view newick) {
    Ancestry cherries;

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
            cherries.push_back({c1, c2, p});

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

    return cherries;
}

Ancestry getCherriesNoParents(std::string_view newick) {
    const size_t newickLength = newick.length();

    Ancestry ancestry;
    ancestry.reserve(newickLength);

    std::vector<int> stack;

    for (size_t i = 0; i < newickLength; ++i) {
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

    // To keep track + update of the smallest descendant of each triplet
    // This allows us to reconstruct the triplets as they should appear
    // using the Phylo2Vec construction
    // Note: the first numLeaves indices are not used
    std::vector<int> minDesc(numNodes, -1);

    // Sort the ancestry by their parent node (ascending order)
    std::qsort(ancestry.data(), numCherries, sizeof(std::array<int, 3>),
               [](const void *a, const void *b) {
                   const auto &arr1 = *static_cast<const std::array<int, 3> *>(a);
                   const auto &arr2 = *static_cast<const std::array<int, 3> *>(b);
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

void orderCherriesNoParents(Ancestry &cherries) {
    std::vector<int> leaves;
    std::unordered_map<int, int> visited;

    for (size_t i = 0; i < cherries.size(); ++i) {
        auto &[c1, c2, cMax] = cherries[i];
        int cMin = std::min(c1, c2);

        int toProcess = cMax;

        if (visited.find(cMin) != visited.end() && visited[cMin] < cMax) {
            toProcess = visited[cMin];
        }

        leaves.push_back(toProcess);
        visited[cMin] = toProcess;
    }

    // argsort with descending order
    // argsort used here for to_matrix to reorder BLs
    std::vector<size_t> indices(cherries.size());
    std::iota(indices.begin(), indices.end(), 0);
    // stable sort is important to keep the order of cherries
    std::stable_sort(indices.begin(), indices.end(),
                     [&leaves](size_t i, size_t j) { return leaves[i] > leaves[j]; });

    // Reorder the cherries in place using the sorted indices
    Ancestry temp;
    temp.reserve(cherries.size());
    for (size_t i : indices) {
        temp.push_back(cherries[i]);
    }
    cherries = std::move(temp);
}

PhyloVec buildVector(Ancestry cherries) {
    const size_t numCherries = cherries.size();
    const size_t numLeaves = numCherries + 1;

    PhyloVec v(numCherries, 0);

    FenwickTree bit(numLeaves);

    // Note: v[0] is always 0
    // but starting with i = 1 makes some tests fail (weird)
    for (size_t i = 0; i < numCherries; ++i) {
        auto &[c1, c2, cMax] = cherries[i];

        unsigned int idx = bit.prefix_sum(cMax - 1);

        // Reminder: v[i] = j --> branch i yields leaf j
        v[cMax - 1] = idx == 0 ? std::min(c1, c2) : cMax - 1 + idx;

        bit.update(cMax, 1);
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
    Ancestry ancestry = getCherriesNoParents(newick.substr(0, newick.length() - 1));

    orderCherriesNoParents(ancestry);

    PhyloVec v = buildVector(ancestry);

    return v;
}