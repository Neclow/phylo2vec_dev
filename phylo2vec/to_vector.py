import re

import numpy as np

POS_INT_PATTERN = re.compile(r"\b\d+\b")


def _remove_branch_length_annotations(newick):
    """Remove branch lengths annotations from a Newick string

    Example: "(((2:0.02,1:0.01),0:0.041),3:1.42);" --> "(((2,1),0),3);"

    Parameters
    ----------
    newick : str
        Newick representation of a tree
    """
    pass


def _remove_parent_annotations(newick):
    """Remove parent nodes from a Newick string

    Example: "(((2,1)4,0)5,3)6;" --> "(((2,1),0),3);"

    Parameters
    ----------
    newick : str
        Newick representation of a tree
    """
    pass


def to_vector(newick, n_leaves=None, verbose=False):
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
        # n_leaves = max(int(s) for s in re.findall(POS_INT_PATTERN, nw)) + 1
        n_leaves = len(re.findall(POS_INT_PATTERN, newick))

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
