"""
Methods to convert a Newick to a Phylo2Vec vector.

Two main methods:
    - to_vector for a Newick with parent labels
    - to_vector_no_parents for a Newick without parent labels
"""

import numba as nb
import numpy as np


def _reduce(newick):
    ancestry = []

    def do_reduce(ancestry, newick):
        for i, char in enumerate(newick):
            if char == "(":
                open_idx = i + 1
            elif char == ")":
                child1, child2 = newick[open_idx:i].split(",", 2)
                parent = newick[i + 1 :].split(",", 1)[0].split(")", 1)[0]

                ancestry.append(
                    [
                        int(child1),
                        int(child2),
                        int(parent),
                    ]
                )
                newick = newick[: open_idx - 1] + newick[i + 1 :]

                return do_reduce(ancestry, newick)

    do_reduce(ancestry, newick[:-1])

    return np.array(ancestry, dtype=np.int16)


def _reduce_no_parents(newick):
    ancestry = []

    def do_reduce(ancestry, newick):
        for i, char in enumerate(newick):
            if char == "(":
                open_idx = i + 1
            elif char == ")":
                child1, child2 = newick[open_idx:i].split(",", 2)

                child1 = int(child1)
                child2 = int(child2)

                ancestry.append([child1, child2, max(child1, child2)])

                newick = newick.replace(
                    newick[open_idx - 1 : i + 1], f"{min(child1, child2)}"
                )

                return do_reduce(ancestry, newick)

    do_reduce(ancestry, newick[:-1])

    return np.array(ancestry, dtype=np.int16)


@nb.njit(cache=True)
def _order_cherries(ancestry):
    ancestry_sorted = ancestry[np.argsort(ancestry[:, -1]), :]

    small_children = nb.typed.Dict.empty(
        key_type=nb.types.int16, value_type=nb.types.int16
    )

    for i, row in enumerate(ancestry_sorted):
        c1, c2, p = row

        parent_c1, parent_c2 = small_children.get(c1, c1), small_children.get(c2, c2)

        small_children[p] = min(parent_c1, parent_c2)

        ancestry_sorted[i, :] = [parent_c1, parent_c2, max(parent_c1, parent_c2)]
    return ancestry_sorted


@nb.njit(cache=True)
def _order_cherries_no_parents(cherries):
    n_cherries = cherries.shape[0]

    old_cherries = cherries.copy()

    idxs = np.zeros((n_cherries,), dtype=np.uint8)

    for i in range(n_cherries):
        unvisited = np.ones((n_cherries + 1,), dtype=np.uint8)
        max_leaf = -1

        for j, ch in enumerate(old_cherries):
            if idxs[j] == 1:
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
        idxs[idx] = 1

    return cherries


@nb.njit(cache=True)
def _build_vector(cherries):
    v_res = np.zeros((cherries.shape[0],), dtype=np.uint16)
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
    """Convert a Newick string with parent labels to a vector

    Parameters
    ----------
    newick : str
        Newick string for a tree

    Returns
    -------
    v : numpy.ndarray
        Phylo2Vec vector
    """
    ancestry = _reduce(newick)

    cherries = _order_cherries(ancestry)

    v = _build_vector(cherries)

    return v


def to_vector_no_parents(newick_no_parents):
    """Convert a Newick string without parent labels to a vector

    Parameters
    ----------
    newick_no_parents : str
        Newick string for a tree

    Returns
    -------
    v : numpy.ndarray
        Phylo2Vec vector
    """
    ancestry = _reduce_no_parents(newick_no_parents)

    cherries = _order_cherries_no_parents(ancestry)

    v = _build_vector(cherries)

    return v
