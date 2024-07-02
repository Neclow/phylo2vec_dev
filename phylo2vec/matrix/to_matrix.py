import numpy as np

from phylo2vec.base.to_vector import (
    _build_vector,
    _find_cherries,
    _order_cherries_no_parents,
)


def _reduce_with_bls(newick):
    ancestry = []
    bls = []

    def do_reduce(ancestry, bls, newick):
        for i, char in enumerate(newick):
            if char == "(":
                open_idx = i + 1
            elif char == ")":
                child1, child2 = newick[open_idx:i].split(",", 2)
                parent = newick[i + 1 :].split(",", 1)[0].split(")", 1)[0]

                child1, bl1 = child1.split(":")
                child2, bl2 = child2.split(":")
                parent = parent.split(":")[0]

                ancestry.append(
                    [
                        int(child1),
                        int(child2),
                        int(parent),
                    ]
                )

                bls.append([float(bl1), float(bl2)])
                newick = newick[: open_idx - 1] + newick[i + 1 :]

                return do_reduce(ancestry, bls, newick)

    do_reduce(ancestry, bls, newick[:-1])

    return np.array(ancestry, dtype=np.int16), np.array(bls, dtype=np.float16)


def _reduce_no_parents_with_bls(newick):
    ancestry = []

    bls = []

    def do_reduce(ancestry, bls, newick):
        for i, char in enumerate(newick):
            if char == "(":
                open_idx = i + 1
            elif char == ")":
                child1, child2 = newick[open_idx:i].split(",", 2)

                child1, bl1 = child1.split(":")
                child2, bl2 = child2.split(":")

                child1 = int(child1)
                child2 = int(child2)

                ancestry.append([child1, child2, max(child1, child2)])

                bls.append([bl1, bl2])

                newick = newick.replace(
                    newick[open_idx - 1 : i + 1], f"{min(child1, child2)}"
                )

                return do_reduce(ancestry, bls, newick)

    do_reduce(ancestry, bls, newick[:-1])

    return np.array(ancestry, dtype=np.int16), np.array(bls, dtype=np.float16)


def to_matrix(newick):
    ancestry, bls = _reduce_with_bls(newick)

    cherries = _find_cherries(ancestry)

    v = _build_vector(cherries)

    return v, bls


def to_matrix_no_parents(newick_no_parents):
    ancestry, bls = _reduce_no_parents_with_bls(newick_no_parents)

    cherries = _order_cherries_no_parents(ancestry)

    v = _build_vector(cherries)

    return v, bls