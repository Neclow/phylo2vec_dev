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
            stack.push_back(p);
        } else if (c >= '0' && c <= '9') {
            // Get the next node and push it the stack
            size_t end;
            int n = stoi_substr(newick, i, &end);
            stack.push_back(n);
            i = end - 1;
        }
    }

    return ancestry;
}

Ancestry getCherriesNoParents(std::string_view newick) {
    Ancestry ancestry;

    std::vector<int> stack;

    for (size_t i = 0; i < newick.length(); ++i) {
        char c = newick[i];

        if (c == ')') {
            int c2 = stack.back();
            stack.pop_back();

            int c1 = stack.back();
            stack.pop_back();

            int cMax = std::max(c1, c2);

            ancestry.push_back({c1, c2, cMax});
            stack.push_back(std::min(c1, c2));
        } else if (c >= '0' && c <= '9') {
            size_t end;
            int node = stoi_substr(newick, i, &end);
            stack.push_back(node);
            i = end - 1;
        }
    }

    return ancestry;
}

void orderCherries(Ancestry &ancestry) {
    const size_t numCherries = ancestry.size();

    std::qsort(ancestry.data(), numCherries, sizeof(std::array<int, 3>),
               [](const void *a, const void *b) {
                   const auto &arr1 =
                       *static_cast<const std::array<int, 3> *>(a);
                   const auto &arr2 =
                       *static_cast<const std::array<int, 3> *>(b);
                   return arr1[2] - arr2[2];
               });

    std::vector<int> child_min(numCherries * 2 + 1, -1);

    for (size_t i = 0; i < numCherries; ++i) {
        auto &[c1, c2, p] = ancestry[i];

        int parentOfChild1 = child_min[c1] != -1 ? child_min[c1] : c1;
        int parentOfChild2 = child_min[c2] != -1 ? child_min[c2] : c2;

        child_min[p] = std::min(parentOfChild1, parentOfChild2);

        ancestry[i] = {parentOfChild1, parentOfChild2,
                       std::max(parentOfChild1, parentOfChild2)};
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
        std::vector<std::uint8_t> d(2 * numCherries, 0);
        int maxLeaf = -1;

        int idx;

        for (int j = 0; j < numCherries; ++j) {
            if (idxs[j]) {
                continue;
            }

            auto &[c1, c2, _] = oldAncestry[j];

            if (!(d[c1] || d[c2])) {
                int cMax = std::max(c1, c2);
                if (cMax > maxLeaf) {
                    maxLeaf = std::max(c1, c2);
                    idx = j;
                }
            }

            d[c1] = 1;
            d[c2] = 1;
        }

        idxs[idx] = 1;

        ancestry[i] = oldAncestry[idx];
    }
}

PhyloVec buildVector(Ancestry const &cherries) {
    const size_t numCherries = cherries.size();

    PhyloVec v(numCherries);

    for (int i = numCherries - 1; i >= 0; --i) {
        auto &[c1, c2, p] = cherries[i];

        int cMin = std::min(c1, c2);
        int cMax = std::max(c1, c2);

        int idx = 0;

        // Find first row containing cMax
        for (size_t j = 0; j < numCherries; j++) {
            if (cherries[j][2] <= cMax) {
                if (cherries[j][0] == cMax || cherries[j][1] == cMax) {
                    break;
                } else {
                    ++idx;
                }
            }
        }

        v[cMax - 1] = idx == 0 ? cMin : cMax - 1 + idx;
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