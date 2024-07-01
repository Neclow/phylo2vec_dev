#include "to_newick.hpp"

#include <iostream>
#include <unordered_map>

Ancestry getAncestry(const PhyloVec &v)
{
    // Matrix with 3 columns: child1, child2, parent
    Ancestry ancestry(v.size());

    // This is the first pair, we start with (0, 1)
    std::vector<std::pair<int, int>> pairs;

    pairs.push_back(std::make_pair(0, 1));

    int next_leaf;

    // The goal here is to add mergers like in the previous iteration
    for (std::size_t i = 1; i < v.size(); ++i)
    {
        next_leaf = i + 1;

        if (v[i] <= i)
        {
            /*
            If v[i] <= i, it's an easy BD
            We now that the next pair to add now is (v[i], next_leaf)
            (as the branch leading to v[i] gives birth to the next_leaf)
            Why pairs.insert(0)? Let's take an example with [0, 0]
            We initially have (0, 1), but 0 gives birth to 2 afterwards
            So the "shallowest" pair is (0, 2)
            */
            pairs.insert(pairs.begin(), std::make_pair(v[i], next_leaf));
        }
        else
        {
            /*
            If v[i] > i, it's not the branch leading v[i] that gives birth but an internal branch
            Remark: it will not be the "shallowest" pair, so we do not insert it at position 0
            len(pairs) = number of pairings we did so far
            So what v[i] - len(pairs) gives us is the depth of the next pairing
            And pairs[v[i] - len(pairs) - 1][0] is a node that we processed beforehand
            which is deeper than the branch v[i]
            */
            pairs.insert(
                pairs.begin() + v[i] - pairs.size(),
                std::make_pair(pairs[v[i] - pairs.size() - 1].first, next_leaf));
        }
    }

    // Dictionary to keep track of the following relationship: child->highest parent
    std::unordered_map<int, int> parents;

    // Dictionary to keep track of siblings (i.e., sister nodes)
    std::unordered_map<int, int> siblings;

    // Leaves are number 0, 1, ..., n_leaves - 1, so the next parent is n_leaves
    int next_parent = v.size() + 1;

    int child1, child2;
    int parent_child1, parent_child2;
    int sibling_child1, sibling_child2;

    for (std::size_t i = 0; i < pairs.size(); ++i)
    {
        child1 = pairs[i].first;
        child2 = pairs[i].second;

        parent_child1 = parents.find(child1) != parents.end() ? parents[child1] : child1;
        parent_child2 = parents.find(child2) != parents.end() ? parents[child2] : child2;

        ancestry[i] = {parent_child1, parent_child2, next_parent};

        // Change the parents of the current children
        parents[child1] = next_parent;
        parents[child2] = next_parent;

        sibling_child1 = siblings.find(child1) != siblings.end() ? siblings[child1] : child1;
        sibling_child2 = siblings.find(child2) != siblings.end() ? siblings[child2] : child2;

        // Change the parents of the siblings
        parents[sibling_child1] = next_parent;
        parents[sibling_child2] = next_parent;

        // Change the previous parents of the child if there are any
        parents[parent_child1] = next_parent;
        parents[parent_child1] = next_parent;

        // Update siblings
        ++next_parent;
    }

    return ancestry;
}

std::string buildNewick(const Ancestry &ancestry)
{
    // Row with 2 children of root + root node
    std::array<int, 3> row = ancestry.back();
    int c1 = row[0];
    int c2 = row[1];
    int p = row[2];

    std::string c1_str = std::to_string(c1);

    std::string newick = "(" + c1_str + "," + std::to_string(c2) + ")" + std::to_string(p) + ";";

    std::unordered_map<int, int> node_idxs;
    node_idxs = {{c1, 1},
                 {c2, 2 + c1_str.length()}};

    std::vector<int> queue;

    size_t n_max = ancestry.size();

    if (c1 > n_max)
    {
        queue.push_back(c1);
    }
    if (c2 > n_max)
    {
        queue.push_back(c2);
    }

    int next_parent;
    std::string sub_newick;
    std::string p_str;

    for (std::size_t i = 1; i < ancestry.size(); ++i)
    {
        next_parent = queue.back();

        queue.pop_back();

        row = ancestry[next_parent - n_max - 1];
        int c1 = row[0];
        int c2 = row[1];
        int p = row[2];

        c1_str = std::to_string(c1);
        p_str = std::to_string(p);

        sub_newick = "(" + c1_str + "," + std::to_string(c2) + ")" + p_str;

        newick = newick.substr(0, node_idxs[p]) + sub_newick + newick.substr(node_idxs[p] + p_str.length());

        node_idxs[c1] = node_idxs[p] + 1;
        node_idxs[c2] = node_idxs[c1] + 1 + c1_str.length();

        if (c1 > n_max)
        {
            queue.push_back(c1);
        }
        if (c2 > n_max)
        {
            queue.push_back(c2);
        }
    }

    return newick;
}

std::string toNewick(const PhyloVec &v)
{
    return buildNewick(getAncestry(v));
}