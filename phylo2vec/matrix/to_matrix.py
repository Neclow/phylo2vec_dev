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
            # Add the pair (bl1, bl2)
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


# TODO: handle newick without parents + parent BLs?
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

            if c1 < c2:
                c_min = c1
                c_max = c2
            else:
                c_min = c2
                c_max = c1

            c_min, c_max = sorted([c1, c2])

            # No parent annotation --> store the max leaf
            ancestry.append([c1, c2, c_max])
            bls.append([bl1, bl2])

            # Push the min leaf to the stack
            annotated_node, end = _node_substr(newick, i)
            stack.append(c_min)
            try:
                bl_stack.append(float(annotated_node[1:]))
                i = end - 1
            except ValueError:
                break

        elif "0" <= char <= "9":
            annotated_node, end = _node_substr(newick, i)

            try:
                # print(annotated_node)
                node, bln = annotated_node.split(":", 1)

                stack.append(int(node))
                bl_stack.append(float(bln))
            except ValueError:
                stack.append(min(stack[-2:]))
                bl_stack.append(float(annotated_node))

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
        # e.g., ((1:0.1,2:0.2):0.6,(3:0.4,4:0.5):0.7);

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
