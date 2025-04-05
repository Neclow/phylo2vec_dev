"""
Methods to convert a Newick to a Phylo2Mat matrix.

Two main methods:
    - to_matrix for a Newick with parent labels + branch lengths
    - to_matrix_no_parents for a Newick without parent labels + branch lengths
"""

import numpy as np

from phylo2vec.base.to_vector import (
    _build_vector,
    _order_cherries,
    _order_cherries_no_parents,
    _node_substr,
)


def _get_cherries_with_bls(newick):
    ancestry = []
    bls = []
    stack = []
    bl_stack = []

    i = 0
    while i < len(newick):
        char = newick[i]
        if char == ")":
            i += 1

            # Pop the children nodes from the stack
            c2 = stack.pop()
            c1 = stack.pop()

            # Pop the BLs from the BL stack
            bl2 = bl_stack.pop()
            bl1 = bl_stack.pop()
            bls.append([bl1, bl2])

            # Get the parent node after )
            annotated_p, end = _node_substr(newick, i)
            i = end - 1

            # Add the triplet (c1, c2, p)
            try:
                p, blp = annotated_p.split(":", 1)
                ancestry.append([c1, c2, int(p)])

                # Push the parent node to the stack
                # Push the parent BL to the BL stack
                stack.append(int(p))
                bl_stack.append(float(blp))
            except ValueError:
                # No BL for the root, so we just add the triplet
                ancestry.append([c1, c2, int(annotated_p)])

        elif "0" <= char <= "9":
            # Get the next node and push it to the stack
            annotated_node, end = _node_substr(newick, i)

            node, bln = annotated_node.split(":", 1)

            stack.append(int(node))
            bl_stack.append(float(bln))

            i = end - 1

        i += 1

    return np.asarray(ancestry, dtype=np.int32), np.asarray(bls, dtype=np.float16)


# TODO: handle newick without parents & without parent BLs?
# e.g., ((1:0.1,2:0.2),(3:0.4,4:0.5));
def _get_cherries_no_parents_with_bls(newick):
    ancestry = []
    bls = []
    stack = []
    bl_stack = []

    i = 0
    while i < len(newick):
        char = newick[i]
        if char == ")":
            i += 1

            # Pop the children nodes from the stack
            c2 = stack.pop()
            c1 = stack.pop()

            # Pop the BLs from the BL stack
            bl2 = bl_stack.pop()
            bl1 = bl_stack.pop()

            c_min, c_max = sorted([c1, c2])

            # No parent annotation --> store the max leaf
            ancestry.append([c1, c2, c_max])
            bls.append([bl1, bl2])

            # Find the parental BL
            # Ex: ":0.2"
            annotated_node, end = _node_substr(newick, i)
            i = end - 1

            if len(annotated_node) == 0 and end == len(newick) - 1:
                # if this is true, we reached the root which we don't have a parent BL
                # We could break, but "continue" might prevent silent errors
                break

            # Push the min leaf to the stack
            stack.append(c_min)
            # Push the parent BL to the BL stack
            bl_stack.append(float(annotated_node[1:]))

        elif "0" <= char <= "9":
            annotated_node, end = _node_substr(newick, i)

            node, bln = annotated_node.split(":", 1)

            stack.append(int(node))
            bl_stack.append(float(bln))

            i = end - 1

        i += 1

    return np.asarray(ancestry, dtype=np.int32), np.asarray(bls, dtype=np.float16)


def to_matrix(newick):
    """
    Convert a Newick string with parent labels and branch lengths to a matrix

    This functions wraps:
      * an adaptation of ```_get_cherries```with branch lengths
      * the base function ```_order_cherries```
      * the base function ```_build_vector```

    Parameters
    ----------
    newick : str
        Newick string for a tree
        # e.g., ((0:0.1,1:0.2)4:0.6,(2:0.4,3:0.5)5:0.7);

    Returns
    -------
    m : numpy.ndarray
        Phylo2Mat matrix
    """
    ancestry, bls = _get_cherries_with_bls(newick)

    cherries, idxs = _order_cherries(ancestry)

    bls = bls[idxs]

    v = _build_vector(cherries)

    m = np.concatenate([v[:, None], bls.astype(np.float16)], axis=1)

    return m


def to_matrix_no_parents(newick_no_parents):
    """
    Convert a Newick string with parent labels and branch lengths to a matrix

    This functions wraps:
      * an adaptation of ```_get_cherries_no_parents```with branch lengths
      * the base function ```_order_cherries_no_parents```
      * the base function ```_build_vector```

    Parameters
    ----------
    newick : str
        Newick string for a tree
        # e.g., ((0:0.1,1:0.2):0.6,(3:0.4,4:0.5):0.7);

    Returns
    -------
    m : numpy.ndarray
        Phylo2Mat matrix
    """
    ancestry, bls = _get_cherries_no_parents_with_bls(newick_no_parents)

    cherries, idxs = _order_cherries_no_parents(ancestry)

    bls = bls[idxs]

    v = _build_vector(cherries)

    m = np.concatenate([v[:, None], bls.astype(np.float16)], axis=1)

    return m
