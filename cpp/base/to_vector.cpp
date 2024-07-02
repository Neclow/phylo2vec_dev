#include "to_vector.hpp"
#include <algorithm>
#include <cmath>
// // #include <unordered_map>

// size_t getNumLeavesFromNewick(std::string const &newick) {
//     size_t idx = newick.find_last_of(")");

//     size_t root = std::stoi(newick.substr(idx + 1, newick.length() - idx));
//     return root / 2 + 1;
// }

void doReduce(Ancestry &ancestry, std::string &newick) {
    size_t openIdx, commaIdx;

    int c1, c2, p;

    std::string children, parent, rest;

    for (size_t i = 0; i < newick.length(); ++i) {
        if (newick[i] == '(') {
            openIdx = i + 1;
        } else if (newick[i] == ')') {
            children = newick.substr(openIdx, i - openIdx);

            commaIdx = children.find(',');

            c1 = std::stoi(children.substr(0, commaIdx));
            c2 = std::stoi(children.substr(commaIdx + 1));

            rest = newick.substr(i + 1);

            parent = rest.substr(0, rest.find_first_of(",)"));

            p = std::stoi(parent);

            ancestry.push_back({c1, c2, p});

            newick = newick.substr(0, openIdx - 1) + rest;

            return doReduce(ancestry, newick);
        }
    }
}

Ancestry reduce(std::string &newick) {
    Ancestry ancestry;

    newick = newick.substr(0, newick.length() - 1);

    doReduce(ancestry, newick);

    return ancestry;
}

void toCherries(Ancestry &ancestry) {
    std::qsort(ancestry.data(), ancestry.size(), sizeof(std::array<int, 3>),
               [](const void *a, const void *b) {
                   const auto &arr1 =
                       *static_cast<const std::array<int, 3> *>(a);
                   const auto &arr2 =
                       *static_cast<const std::array<int, 3> *>(b);
                   return arr1[2] - arr2[2];
               });

    std::vector<int> child_min(ancestry.size() * 2 + 1, -1);

    int parentOfChild1;
    int parentOfChild2;

    for (size_t i = 0; i < ancestry.size(); ++i) {
        auto &[c1, c2, p] = ancestry[i];

        parentOfChild1 = child_min[c1] != -1 ? child_min[c1] : c1;
        parentOfChild2 = child_min[c2] != -1 ? child_min[c2] : c2;

        // parentofChild1 = child_min.find(c1) != child_min.end() ?
        // child_min[c1] : c1;
        // parentOfChild2 = child_min.find(c2) != child_min.end() ?
        // child_min[c2] : c2;

        child_min[p] = std::min(parentOfChild1, parentOfChild2);

        ancestry[i] = {parentOfChild1, parentOfChild2,
                       std::max(parentOfChild1, parentOfChild2)};
    }
}

PhyloVec buildVector(Ancestry const &cherries) {
    size_t numLeaves = cherries.size();

    PhyloVec v(numLeaves);

    int c_max;

    int idx;

    std::vector<std::array<int, 2>> subset;
    subset.reserve(cherries.size());

    for (int i = numLeaves - 1; i >= 0; --i) {
        auto &[c1, c2, p] = cherries[i];

        c_max = std::max(c1, c2);

        for (size_t j = 0; j < numLeaves; j++) {
            if (cherries[j][2] <= c_max) {
                subset.push_back({cherries[j][0], cherries[j][1]});
            }
        }

        for (size_t j = 0; j < subset.size(); j++) {
            if (subset[j][0] == c_max || subset[j][1] == c_max) {
                idx = static_cast<int>(j);
                break;
            }
        }

        v[c_max - 1] = idx == 0 ? std::min(c1, c2) : c_max - 1 + idx;

        subset.clear();
    }

    return v;
}

PhyloVec toVector(std::string &newick) {
    Ancestry ancestry = reduce(newick);

    toCherries(ancestry);

    PhyloVec v = buildVector(ancestry);

    return v;
}