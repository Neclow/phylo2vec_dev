import numba as nb

from phylo2vec.base.to_newick import _get_ancestry


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


@nb.njit(cache=True)
def to_newick(v, bls):
    ancestry = _get_ancestry(v)

    newick = _build_newick_with_bls(ancestry, bls)

    return newick
