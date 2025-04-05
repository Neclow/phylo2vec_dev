"""
Methods to convert a Newick to a Phylo2Vec vector.

Two main methods:
    - to_vector for a Newick with parent labels
    - to_vector_no_parents for a Newick without parent labels
"""

import numba as nb
import numpy as np


def _node_substr(s, start):
    end = 0

    for i, c in enumerate(s[start:], start):
        if c in (",", ")", ";"):
            end = i
            break

    if end == 0:
        raise ValueError(f"Failed to parse node in Newick string: {s}")

    s_substr = s[start:end]

    return s_substr, end


def _get_cherries(newick):
    ancestry = []
    stack = []

    i = 0
    while i < len(newick):
        char = newick[i]
        if char == ")":
            i += 1

            # Pop the children nodes from the stack
            c2 = stack.pop()
            c1 = stack.pop()

            # Get the parent node after )
            p, end = _node_substr(newick, i)
            i = end - 1

            # Add the triplet (c1, c2, p)
            ancestry.append([c1, c2, int(p)])

            # Push the parent node to the stack
            stack.append(int(p))
        elif "0" <= char <= "9":
            # Get the next node and push it to the stack
            node, end = _node_substr(newick, i)

            stack.append(int(node))

            i = end - 1

        i += 1

    return np.asarray(ancestry, dtype=np.int32)


def _get_cherries_no_parents(newick):
    ancestry = []
    stack = []

    i = 0
    while i < len(newick):
        char = newick[i]
        if char == ")":
            # Pop the children nodes from the stack
            c2 = stack.pop()
            c1 = stack.pop()

            c_min, c_max = sorted([c1, c2])

            # No parent annotation --> store the max leaf
            ancestry.append([c1, c2, c_max])

            # Push the min leaf to the stack
            stack.append(c_min)
        elif "0" <= char <= "9":
            # Get the next leaf and push it to the stack
            node, end = _node_substr(newick, i)

            stack.append(node)

            i = end - 1

        i += 1

    return np.asarray(ancestry, dtype=np.int32)


@nb.njit(cache=True)
def _order_cherries(ancestry):
    idxs = np.argsort(ancestry[:, -1])

    ancestry_sorted = ancestry[idxs, :]

    small_children = nb.typed.Dict.empty(
        key_type=nb.types.int32, value_type=nb.types.int32
    )

    for i, row in enumerate(ancestry_sorted):
        c1, c2, p = row

        parent_c1, parent_c2 = small_children.get(c1, c1), small_children.get(c2, c2)

        small_children[p] = min(parent_c1, parent_c2)

        ancestry_sorted[i, :] = [parent_c1, parent_c2, max(parent_c1, parent_c2)]

    return ancestry_sorted, idxs


@nb.njit(cache=True)
def _order_cherries_no_parents(cherries):
    n_cherries = cherries.shape[0]

    old_cherries = cherries.copy()

    idxs = []  # np.zeros((n_cherries,), dtype=np.uint8)

    for i in range(n_cherries):
        unvisited = np.ones((n_cherries + 1,), dtype=np.uint8)
        max_leaf = -1
        idx = -1

        for j, ch in enumerate(old_cherries):
            if j in idxs:
                continue

            c1, c2, c_max = ch

            if unvisited[c1] and unvisited[c2]:
                if c_max > max_leaf:
                    max_leaf = c_max
                    idx = j

            unvisited[c1] = 0
            unvisited[c2] = 0

        # Swap the rows for the new ancestry
        # row idx becomes row i
        cherries[i] = old_cherries[idx]

        # Row idx has been processed
        # idxs[idx] = 1
        idxs.append(idx)

    return cherries, np.asarray(idxs)


@nb.njit(cache=True)
def _build_vector(cherries):
    v_res = np.zeros((cherries.shape[0],), dtype=np.uint32)
    for i in range(cherries.shape[0]):
        c1, c2, c_max = cherries[i]

        idx = 0

        for j in range(i):
            if cherries[j][-1] <= c_max:
                idx += 1

        if idx == 0:
            v_res[c_max - 1] = min(c1, c2)
        else:
            v_res[c_max - 1] = c_max - 1 + idx

    return v_res


def to_vector(newick):
    """
    Convert a Newick string with parent labels to a vector

    Parameters
    ----------
    newick : str
        Newick string for a tree

    Returns
    -------
    v : numpy.ndarray
        Phylo2Vec vector
    """
    ancestry = _get_cherries(newick)

    cherries, _ = _order_cherries(ancestry)

    v = _build_vector(cherries)

    return v


def to_vector_no_parents(newick_no_parents):
    """
    Convert a Newick string without parent labels to a vector

    Parameters
    ----------
    newick_no_parents : str
        Newick string for a tree

    Returns
    -------
    v : numpy.ndarray
        Phylo2Vec vector
    """
    ancestry = _get_cherries_no_parents(newick_no_parents)

    cherries, _ = _order_cherries_no_parents(ancestry)

    v = _build_vector(cherries)

    return v
