import warnings

import numba as nb
import numpy as np


@nb.njit(cache=True)
def _get_ancestry(v):
    """

    New get ancestry, but with a MUCH clearer interpretation of what v[i] is.
    v[i] = which BRANCH we do the pairing from


    The initial situation looks like this:
                      R
                      |
                      | --> branch 2
                    // \\
      branch 0 <-- //   \\  --> branch 1
                   0     1

    For v[1], we have 3 possible branches too choose from.
    v[1] = 0 or 1 indicates that we branch out from branch 0 or 1, respectively.
    The new branch yields leaf 2 (like in ordered trees)

    v[1] = 2 is somewhat similar: we create a new branch from R that yields leaf 2

    Parameters
    ----------
    v : list
        Phylo2Vec

    Returns
    -------
    ancestry : numpy.array
        Ancestry
    """
    # This is the first pair, we start with (0, 1)
    pairs = [(0, 1)]

    # The goal here is to add mergers like in the previous iteration
    for i in range(1, len(v)):
        next_leaf = i + 1

        if v[i] <= i:
            # If v[i] <= i, it's an easy BD
            # We now that the next pair to add now is (v[i], next_leaf)
            # (as the branch leading to v[i] gives birth to the next_leaf)
            # Why pairs.insert(0)? Let's take an example with [0, 0]
            # We initially have (0, 1), but 0 gives birth to 2 afterwards
            # So the "deepest" pair is (0, 2)
            pairs.insert(0, (v[i], next_leaf))
        else:
            # If v[i] > i, it's not the branch leading v[i] that gives birth but an internal branch
            # Remark 1: it will not be the "deepest" pair, so we do not insert it at position 0
            # len(pairs) = number of pairings we did so far
            # So what v[i] - len(pairs) gives us is the depth of the next pairing
            # And pairs[v[i] - len(pairs) - 1][0] is a node that we processed beforehand
            # which is deeper than the branch v[i]
            pairs.insert(
                v[i] - len(pairs), (pairs[v[i] - len(pairs) - 1][0], next_leaf)
            )

    # We have our pairs, we can now build our ancestry
    # Matrix with 3 columns: child1, child2, parent
    ancestry = np.zeros((len(pairs), 3), dtype=np.int32)

    # Dictionary to keep track of the following relationship: child->highest parent
    parents = nb.typed.Dict.empty(key_type=nb.types.int64, value_type=nb.types.int64)

    # Dictionary to keep track of siblings (i.e., sister nodes)
    siblings = nb.typed.Dict.empty(key_type=nb.types.int64, value_type=nb.types.int64)

    # Leaves are number 0, 1, ..., n_leaves - 1, so the next parent is n_leaves
    next_parent = len(v) + 1

    for i, pair in enumerate(pairs):
        child1, child2 = pair

        parent_child1 = parents.get(child1, child1)
        parent_child2 = parents.get(child2, child2)

        ancestry[i, :] = [parent_child1, parent_child2, next_parent]

        # Change the parents of the current children
        parents[child1] = next_parent
        parents[child2] = next_parent

        # Change the parents of the siblings
        parents[siblings.get(child1, child1)] = next_parent
        parents[siblings.get(child2, child2)] = next_parent

        # Change the previous parents of the child if there are any
        parents[parent_child1] = next_parent
        parents[parent_child2] = next_parent

        # Update siblings
        siblings[child1] = child2
        siblings[child2] = child1

        next_parent += 1

    # ancestry = np.flip(ancestry)

    return ancestry


def _build_newick_old(ancestry):
    """Build a tree from an "ancestry" array

    The input M should always be 3-dimensional with the following format:
    1st column: child 1 parent node
    2nd column: child 2
    3rd column: parent node

    M is processed such that we iteratively write a Newick string
    to describe the tree.

    Parameters
    ----------
    M : numpy.ndarray
        "Ancestry" array of size (n_leaves - 1, 3)

    Returns
    -------
    str
        Newick string
    """
    warnings.warn(
        "This function is deprecated and will be removed soon.",
        DeprecationWarning,
        stacklevel=2,
    )

    # List of parent nodes
    parent_nodes = []

    # List of sub-Newicks
    sub_newicks = []

    for i in range(len(ancestry)):
        child1, child2, parent = ancestry[i, :]
        # Case 1: Both children are parent nodes, so we have sub-newicks for them
        if child1 in parent_nodes and child2 in parent_nodes:
            # Find their indices
            idx1 = parent_nodes.index(child1)
            idx2 = parent_nodes.index(child2)

            # Merge the sub-newicks and add the parent node
            sub_newicks[idx1] = f"({sub_newicks[idx1]},{sub_newicks[idx2]}){parent}"

            # Update the parent node for the 1st children
            parent_nodes[idx1] = parent

            # Discard info on 2nd children as merged with the 1st children
            sub_newicks.remove(sub_newicks[idx2])
            parent_nodes.remove(parent_nodes[idx2])

        # Case 2: only the first child is a parent node
        elif child1 in parent_nodes:
            # Find its index
            idx = parent_nodes.index(child1)

            # Update its sub-Newick:
            # Before: (sub_child1.1, sub_child1.2)child_1
            # After: ((sub_child1.1, sub_child1.2)child_1, child_2)parent
            sub_newicks[idx] = "(" + sub_newicks[idx].replace(
                f"{child1}", f"{child1},{child2}){parent}"
            )

            # Update the parent node (first child is now just an internal node)
            parent_nodes[idx] = parent

        # Case 3: only the second child is a parent node (similar to Case 2)
        elif child2 in parent_nodes:
            idx = parent_nodes.index(child2)
            # Before: (sub_child2.1, sub_child2.2)child_2
            # After: ((sub_child2.1, sub_child2.2)child_2, child_1)parent
            sub_newicks[idx] = "(" + sub_newicks[idx].replace(
                f"{child2}", f"{child2},{child1}){parent}"
            )
            parent_nodes[idx] = parent

        # Case 4: the children nodes have not been added yet
        else:
            # Add a new sub-Newick for this triplet
            sub_newicks.append(f"({child1},{child2}){parent}")

            # Append the parent node
            parent_nodes.append(parent)

    # If everything went well, only 1 "sub-newick" should be left, with only 1 parent: the root node
    newick = sub_newicks[0] + ";"

    return newick


@nb.njit(cache=True)
def _build_newick(ancestry):
    """Build a tree from an "ancestry" array

    The input M should always be 3-dimensional with the following format:
    1st column: child 1 parent node
    2nd column: child 2
    3rd column: parent node

    M is processed such that we iteratively write a Newick string
    to describe the tree.

    Parameters
    ----------
    M : numpy.ndarray
        "Ancestry" array of size (n_leaves - 1, 3)

    Returns
    -------
    str
        Newick string
    """
    c1, c2, p = ancestry[-1, :]

    newick = f"({c1},{c2}){p};"

    node_idxs = {c1: 1, c2: 2 + len(f"{c1}")}

    queue = []

    n_leaves = ancestry.shape[0]

    if c1 > n_leaves:
        queue.append(c1)
    if c2 > n_leaves:
        queue.append(c2)

    for _ in range(1, ancestry.shape[0]):
        next_parent = queue.pop()

        c1, c2, p = ancestry[next_parent - n_leaves - 1, :]

        sub_newick = f"({c1},{c2}){p}"

        newick = (
            newick[: node_idxs[p]] + sub_newick + newick[node_idxs[p] + len(f"{p}") :]
        )

        node_idxs[c1] = node_idxs[p] + 1
        node_idxs[c2] = node_idxs[c1] + 1 + len(f"{c1}")

        if c1 > n_leaves:
            queue.append(c1)
        if c2 > n_leaves:
            queue.append(c2)

    return newick


@nb.njit(cache=True)
def to_newick(v):
    ancestry = _get_ancestry(v)

    newick = _build_newick(ancestry)

    return newick
