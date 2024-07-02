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

void doReduce(Ancestry &ancestry, std::string_view newick) {
    std::vector<int> stack;
    for (size_t i = 0; i < newick.length(); i++) {
        char c = newick[i];
        if (c == ')') {
            i++;
            size_t end;
            int n = stoi_substr(newick, i, &end);
            i = end - 1;

            int c2 = stack.back();
            stack.pop_back();
            int c1 = stack.back();
            stack.pop_back();
            ancestry.push_back({c1, c2, n});
            stack.push_back(n);
        } else if (c >= '0' && c <= '9') {
            size_t end;
            int n = stoi_substr(newick, i, &end);
            stack.push_back(n);
            i = end - 1;
        }
    }
}

Ancestry reduce(std::string_view newick) {
    Ancestry ancestry;

    doReduce(ancestry, newick.substr(0, newick.length() - 1));

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
    Ancestry ancestry = reduce(newick);

    toCherries(ancestry);

    PhyloVec v = buildVector(ancestry);

    return v;
}