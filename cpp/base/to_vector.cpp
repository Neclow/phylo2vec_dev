#include "to_vector.hpp"
#include <algorithm>
#include <charconv>
#include <set>
#include <stdexcept>

#include <iostream>

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

Ancestry reduce(std::string_view newick) {
    Ancestry ancestry;

    std::vector<int> stack;

    for (size_t i = 0; i < newick.length(); i++) {
        char c = newick[i];
        if (c == ')') {
            i++;
            size_t end;
            int p = stoi_substr(newick, i, &end);
            i = end - 1;

            int c2 = stack.back();
            stack.pop_back();

            int c1 = stack.back();
            stack.pop_back();

            ancestry.push_back({c1, c2, p});
            stack.push_back(p);
        } else if (c >= '0' && c <= '9') {
            size_t end;
            int p = stoi_substr(newick, i, &end);
            stack.push_back(p);
            i = end - 1;
        }

        // std::cout << i << std::endl;

        // for (auto elm : stack) {
        //     std::cout << elm << " ";
        // }

        // std::cout << std::endl;
    }

    return ancestry;
}

Ancestry reduceNoParents(std::string_view newick) {
    Ancestry ancestry;

    std::vector<int> stack;

    for (size_t i = 0; i < newick.length(); i++) {
        char c = newick[i];

        if (c == ')') {
            // i++;
            // size_t end;
            // int p = stoi_substr(newick, i, &end);
            // i = end - 1;

            int c2 = stack.back();
            stack.pop_back();

            int c1 = stack.back();
            stack.pop_back();

            int p = std::max(c1, c2);

            ancestry.push_back({c1, c2, p});
            stack.push_back(std::min(c1, c2));
        } else if (c >= '0' && c <= '9') {
            size_t end;
            int p = stoi_substr(newick, i, &end);
            stack.push_back(p);
            i = end - 1;
        }

        for (auto elm : stack) {
            std::cout << elm << " ";
        }

        std::cout << std::endl;
    }

    return ancestry;
}

void toCherries(Ancestry &ancestry) {
    size_t numLeaves = ancestry.size();

    std::qsort(ancestry.data(), numLeaves, sizeof(std::array<int, 3>),
               [](const void *a, const void *b) {
                   const auto &arr1 =
                       *static_cast<const std::array<int, 3> *>(a);
                   const auto &arr2 =
                       *static_cast<const std::array<int, 3> *>(b);
                   return arr1[2] - arr2[2];
               });

    std::vector<int> child_min(numLeaves * 2 + 1, -1);

    int parentOfChild1;
    int parentOfChild2;

    for (size_t i = 0; i < numLeaves; ++i) {
        auto &[c1, c2, p] = ancestry[i];

        parentOfChild1 = child_min[c1] != -1 ? child_min[c1] : c1;
        parentOfChild2 = child_min[c2] != -1 ? child_min[c2] : c2;

        child_min[p] = std::min(parentOfChild1, parentOfChild2);

        ancestry[i] = {parentOfChild1, parentOfChild2,
                       std::max(parentOfChild1, parentOfChild2)};
    }
}

void toCherriesNoParents(Ancestry &ancestry) {
    int numCherries = ancestry.size();

    std::vector<int> idxs(numCherries);

    for (int i = 0; i < numCherries; ++i) {
        int maxLeaf = -1;

        std::set<int> d;

        int idx;

        for (int j = 0; j < numCherries; ++j) {
            int c1 = ancestry[j][0];
            int c2 = ancestry[j][1];

            if (std::find(idxs.begin(), idxs.end(), j) != idxs.end()) {
                continue;
            }

            if (!(d.find(c1) != d.end() || d.find(c2) != d.end())) {
                bool c1_is_leaf = c1 <= numCherries;
                bool c2_is_leaf = c2 <= numCherries;
                if (c1_is_leaf && !(c2_is_leaf)) {
                    if (c1 > maxLeaf) {
                        maxLeaf = c1;
                        idx = j;
                    }
                } else if (c2_is_leaf && !(c1_is_leaf)) {
                    if (c2 > maxLeaf) {
                        maxLeaf = c2;
                        idx = j;
                    }
                } else if (c1_is_leaf && c2_is_leaf) {
                    int cMax = std::max(c1, c2);
                    if (cMax > maxLeaf) {
                        maxLeaf = cMax;
                        idx = j;
                    }
                } else {
                    idx = j;
                }
            }

            d.insert(c1);
            d.insert(c2);
        }

        idxs[i] = idx;
    }

    // Apply order determined by idxs
    for (int i = 0; i < numCherries; ++i) {
        ancestry[i] = ancestry[idxs[i]];
    }
}

PhyloVec buildVector(Ancestry const &cherries) {
    size_t numLeaves = cherries.size();

    PhyloVec v(numLeaves);

    int c_max;

    int idx;

    for (int i = numLeaves - 1; i >= 0; --i) {
        auto &[c1, c2, p] = cherries[i];

        c_max = std::max(c1, c2);

        idx = 0;

        for (size_t j = 0; j < numLeaves; j++) {
            if (cherries[j][2] <= c_max) {
                if (cherries[j][0] == c_max || cherries[j][1] == c_max) {
                    break;
                } else {
                    idx++;
                }
            }
        }

        v[c_max - 1] = idx == 0 ? std::min(c1, c2) : c_max - 1 + idx;
    }

    return v;
}

PhyloVec toVector(std::string_view newick) {
    Ancestry ancestry = reduce(newick.substr(0, newick.length() - 1));

    toCherries(ancestry);

    PhyloVec v = buildVector(ancestry);

    return v;
}

PhyloVec toVectorNoParents(std::string_view newick) {
    Ancestry ancestry = reduceNoParents(newick.substr(0, newick.length() - 1));

    toCherriesNoParents(ancestry);

    PhyloVec v = buildVector(ancestry);

    return v;
}