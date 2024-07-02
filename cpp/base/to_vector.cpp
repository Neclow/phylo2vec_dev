#include "to_vector.hpp"

#include <algorithm>
#include <unordered_map>

void doReduce(Ancestry &ancestry, std::string &newick)
{
    size_t open_idx, comma_idx;

    int child1, child2, p;

    std::string children, parent, rest;

    for (size_t i = 0; i < newick.length(); ++i)
    {
        if (newick[i] == '(')
        {
            open_idx = i + 1;
        }
        else if (newick[i] == ')')
        {
            children = newick.substr(open_idx, i - open_idx);

            comma_idx = children.find(',');

            child1 = std::stoi(children.substr(0, comma_idx));
            child2 = std::stoi(children.substr(comma_idx + 1));

            rest = newick.substr(i + 1);

            parent = rest.substr(0, rest.find_first_of(",)"));

            p = std::stoi(parent);

            ancestry.push_back({child1, child2, p});

            newick = newick.substr(0, open_idx - 1) + rest;

            return doReduce(ancestry, newick);
        }
    }
}

Ancestry reduce(std::string &newick)
{
    Ancestry ancestry;

    newick = newick.substr(0, newick.length() - 1);

    doReduce(ancestry, newick);

    return ancestry;
}

void toCherries(Ancestry &ancestry)
{
    std::qsort(ancestry.data(), ancestry.size(), sizeof(std::array<int, 3>), [](const void *a, const void *b)
               {
              const auto& arr1 = *static_cast<const std::array<int, 3>*>(a);
              const auto& arr2 = *static_cast<const std::array<int, 3>*>(b);
              return arr1[2] - arr2[2]; });
    // std::sort(ancestry.begin(), ancestry.end(), [](const auto &a, const auto &b)
    //           { return a[2] < b[2]; });

    std::unordered_map<int, int> child_min;

    int parent_c1;
    int parent_c2;

    for (size_t i = 0; i < ancestry.size(); ++i)
    {
        auto &[c1, c2, p] = ancestry[i];
        parent_c1 = child_min.find(c1) != child_min.end() ? child_min[c1] : c1;
        parent_c2 = child_min.find(c2) != child_min.end() ? child_min[c2] : c2;

        child_min[p] = std::min(parent_c1, parent_c2);

        ancestry[i] = {parent_c1,
                       parent_c2,
                       std::max(parent_c1, parent_c2)};
    }
}

PhyloVec buildVector(const Ancestry &cherries)
{
    PhyloVec v(cherries.size());

    unsigned int c_max;

    unsigned int idx;

    std::vector<std::array<int, 2>> subset;

    for (int i = cherries.size() - 1; i >= 0; --i)
    {
        auto &[c1, c2, p] = cherries[i];

        c_max = std::max(c1, c2);

        for (size_t j = 0; j < cherries.size(); j++)
        {
            if (cherries[j][2] <= c_max)
            {
                subset.push_back({cherries[j][0], cherries[j][1]});
            }
        }

        for (size_t j = 0; j < subset.size(); j++)
        {
            if (subset[j][0] == c_max || subset[j][1] == c_max)
            {
                idx = static_cast<unsigned int>(j);
                break;
            }
        }

        v[c_max - 1] = idx == 0 ? std::min(c1, c2) : c_max - 1 + idx;

        subset.clear();
    }

    return v;
}

PhyloVec toVector(std::string &newick)
{
    Ancestry ancestry = reduce(newick);

    toCherries(ancestry);

    PhyloVec v = buildVector(ancestry);

    return v;
}