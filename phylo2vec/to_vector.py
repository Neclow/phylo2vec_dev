import numba as nb
import numpy as np

from phylo2vec.utils import find_num_leaves


def _remove_branch_length_annotations(newick):
    """Remove branch lengths annotations from a Newick string

    Example: "(((2:0.02,1:0.01),0:0.041),3:1.42);" --> "(((2,1),0),3);"

    Parameters
    ----------
    newick : str
        Newick representation of a tree
    """
    pass


def integerize(newick):
    pass


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

    return np.array(ancestry)


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

    return np.array(ancestry)


@nb.njit
def _find_cherries(ancestry):
    ancestry_sorted = ancestry[np.argsort(ancestry[:, -1]), :]

    small_children = nb.typed.Dict.empty(
        key_type=nb.types.int64, value_type=nb.types.int64
    )

    for i, row in enumerate(ancestry_sorted):
        c1, c2, p = row

        parent_c1, parent_c2 = small_children.get(c1, c1), small_children.get(c2, c2)

        small_children[p] = min(parent_c1, parent_c2)

        ancestry_sorted[i, :] = [parent_c1, parent_c2, max(parent_c1, parent_c2)]
    return ancestry_sorted


@nb.njit
def _order_cherries_no_parents(cherries):
    cherries_ = np.zeros((len(cherries), 4), dtype=np.int16) - 1
    cherries_[:, :-1] = cherries
    cherries_copy = cherries[:, :-1]

    n_leaves = cherries.shape[0]

    for i in range(n_leaves):
        next_internal = len(cherries) + i + 1

        max_leaf = -1

        d = set()

        for j, ch in enumerate(cherries_copy):
            c1, c2 = ch

            if c1 == -1 and c2 == -1:
                continue

            if not (c1 in d or c2 in d):
                leaf_tmp = ch[ch <= n_leaves]

                if leaf_tmp.shape[0] > 0:
                    max_leaf_tmp = leaf_tmp.max()
                    if max_leaf_tmp > max_leaf:
                        max_leaf = max_leaf_tmp
                        idx = j
                else:
                    idx = j
            d.add(c1)
            d.add(c2)

        cherries_[idx, -1] = next_internal

        cherries_copy[idx] = -1

    cherries_ = cherries_[cherries_[:, -1].argsort()]

    return cherries_[:, :-1]


@nb.njit
def _build_vector(cherries):
    v_res = np.zeros((cherries.shape[0],), dtype=np.int16)
    for i in range(cherries.shape[0] - 1, -1, -1):
        c1, c2, _ = cherries[i]

        c_max = max(c1, c2)

        subset = cherries[cherries[:, -1] <= c_max][:, :-1]

        idx = np.where(subset == c_max)[0][0]

        if idx == 0:
            v_res[c_max - 1] = min(c1, c2)
        else:
            v_res[c_max - 1] = c_max - 1 + idx

    return v_res


def to_vector(newick):
    ancestry = _reduce(newick)

    cherries = _find_cherries(ancestry)

    v = _build_vector(cherries)

    return v


def to_vector_no_parents(newick_no_parents):
    ancestry = _reduce_no_parents(newick_no_parents)

    cherries = _order_cherries_no_parents(ancestry)

    v = _build_vector(cherries)

    return v


def to_vector_old(newick, n_leaves=None):
    """Convert a newick-format tree to its v representation

    Parameters
    ----------
    newick : str
        Newick representation of a tree
    n_leaves : int, optional
        Number of leaves, by default None
        (Saves some computation time if given in advance)
    verbose : bool, optional
        If True, print intermediate results

    Returns
    -------
    v: numpy.ndarray
        v representation of nw
    """
    if n_leaves is None:
        # INFO: fastest when n_leaves is used, this is just in case
        # TODO: n_leaves as non-optional argument?
        n_leaves = find_num_leaves(newick)

    # Phylo2Vec vector
    v = np.zeros(n_leaves, dtype=np.int16)

    # Whether each leaf node has been processed or not
    processed = np.zeros(n_leaves, dtype=bool)

    # TODO: documentation
    vmin = np.zeros(n_leaves, dtype=np.int16)
    labels = np.arange(n_leaves, dtype=np.int16)

    try:
        for _ in range(n_leaves - 1):
            # Name of left leaf
            left_leaf = ""

            for i in range(n_leaves):
                if processed[n_leaves - i - 1] == 0:
                    # Find whether the node with the current label has a sister node
                    label = labels[n_leaves - i - 1]

                    # Is label on the left of a newick pair?
                    if newick.find(f"({label},") > -1:
                        left_sep = f"({label},"
                        right_sep = ")"
                        # Sister node = substring between last left_sep and first right_sep
                        left_leaf = newick.rpartition(left_sep)[2].partition(right_sep)[
                            0
                        ]

                    # Is label on the right of a newick pair?
                    elif newick.find(f",{label})") > -1:
                        left_sep = "("
                        right_sep = f",{label})"
                        # Sister node = substring between last left_sep and first right_sep
                        left_leaf = newick.partition(right_sep)[0].rpartition(left_sep)[
                            2
                        ]

                    # Otherwise --> it has no sister node No sister node --> we can skip it
                    else:
                        continue

                    # If the sister substring is an actual digit, we can stop
                    if left_leaf.isdigit():
                        break

                    # Reset the left_leaf if it wasn't a digit
                    left_leaf = ""

            left_leaf_ind = np.arange(len(labels))[labels == int(left_leaf)][0]
            right_leaf = n_leaves - i - 1

            for n in range(right_leaf + 1, n_leaves):
                if not processed[n]:
                    if vmin[n] == 0:
                        vmin[n] = n
                    else:
                        vmin[n] += 1

            labels[left_leaf_ind] = labels.max() + 1

            if vmin[right_leaf] == 0:
                v[right_leaf] = left_leaf_ind
            else:
                v[right_leaf] = vmin[right_leaf]

            # Update the processed vector
            processed[right_leaf] = True

            # Update the Newick string
            newick = newick.replace(
                f"({left_leaf},{labels[right_leaf]})", str(labels[int(left_leaf_ind)])
            )
            newick = newick.replace(
                f"({labels[right_leaf]},{left_leaf})", str(labels[int(left_leaf_ind)])
            )
    except IndexError as e:
        raise IndexError(
            "Have you tried reroot=True? "
            "Are the Newick nodes integers (and not taxa)? "
            "If the error still persists, your tree might be unrooted or non-binary."
        ) from e
    return v
