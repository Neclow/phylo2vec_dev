"""
Methods to convert a Phylo2Mat matrix to a Newick-format string.
"""

import numba as nb
import numpy as np

from phylo2vec.base.to_newick import _get_ancestry
from phylo2vec.utils.validation import check_m


@nb.njit(cache=True)
def _build_newick_with_bls(anc, bls):
    c1, c2, p = anc[-1, :]

    b1, b2 = bls[-1, :]

    newick = f"({c1}:{b1},{c2}:{b2}){p};"

    node_idxs = {c1: 1, c2: 2 + len(f"{c1}:{b1}")}

    queue = []

    n_max = anc.shape[0]

    if c1 > n_max:
        queue.append(c1)
    if c2 > n_max:
        queue.append(c2)

    for _ in range(1, anc.shape[0]):
        next_parent = queue.pop()

        c1, c2, p = anc[next_parent - n_max - 1, :]

        b1, b2 = bls[next_parent - n_max - 1, :]

        sub_newick = f"({c1}:{b1},{c2}:{b2}){p}"
        newick = (
            newick[: node_idxs[p]] + sub_newick + newick[node_idxs[p] + len(f"{p}") :]
        )

        node_idxs[c1] = node_idxs[p] + 1
        node_idxs[c2] = node_idxs[c1] + 1 + len(f"{c1}:{b1}")

        if c1 > n_max:
            queue.append(c1)
        if c2 > n_max:
            queue.append(c2)

    return newick


def to_newick(m):
    """
    Recover a rooted tree (in Newick format) from a Phylo2Mat m

    The functions wraps the base function ```_get_ancestry```
    and a new _build_newick_with_bls which takes into account branch lengths

    Parameters
    ----------
    m : numpy.array
        Phylo2Mat matrix

    * 1st column is v[i]
    * 2nd column is where leaf i branched out from branch v[i]
    * 3rd column is the branch length leading to leaf i

    Returns
    -------
    newick : str
        Newick tree
    """
    check_m(m)

    v, bls = np.split(m, (1,), axis=1)

    ancestry = _get_ancestry(v.astype(np.uint16)[:, 0])

    newick = _build_newick_with_bls(ancestry, bls.round(6).astype(str))

    return newick


@nb.njit(cache=True)
def to_newick2(v, bls):
    """
    Recover a rooted tree (in Newick format) from a Phylo2Vec vector with branch lengths

    Parameters
    ----------
    v : numpy.array
        Phylo2Vec vector
    bls : numpy.array
        Array of branch lengths

    * 1st column of `bls` is where leaf i branched out from branch v[i]
    * 2nd column of `bls` is the branch length leading to leaf i

    Returns
    -------
    newick : str
        Newick tree
    """
    ancestry = _get_ancestry(v)

    newick = _build_newick_with_bls(ancestry, bls)

    return newick
