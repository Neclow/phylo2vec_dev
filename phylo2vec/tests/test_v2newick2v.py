import numpy as np
import pytest

from ete3 import Tree

from phylo2vec.tests.config import MIN_N_LEAVES, MAX_N_LEAVES, N_REPEATS
from phylo2vec.base.to_newick import to_newick
from phylo2vec.base.to_vector import (
    _find_cherries,
    _order_cherries_no_parents,
    _reduce,
    _reduce_no_parents,
    to_vector,
)
from phylo2vec.utils import sample


MIN_N_LEAVES = 5
MAX_N_LEAVES = 200
N_REPEATS = 10


@pytest.mark.parametrize("n_leaves", range(MIN_N_LEAVES, MAX_N_LEAVES + 1))
def test_v2newick2v(n_leaves):
    for _ in range(N_REPEATS):
        v = sample(n_leaves)
        newick = to_newick(v)
        v2 = to_vector(newick)
        assert np.all(v == v2)


@pytest.mark.parametrize("n_leaves", range(MIN_N_LEAVES, MAX_N_LEAVES + 1))
def test_cherries_no_parents(n_leaves):
    for _ in range(N_REPEATS):
        v = sample(n_leaves)
        newick = to_newick(v)
        newick_no_parents = Tree(newick).write(format=9)
        cherries = _find_cherries(_reduce(newick))
        cherries_no_parents = _order_cherries_no_parents(
            _reduce_no_parents(newick_no_parents)
        )

        assert np.array_equal(cherries, cherries_no_parents)


if __name__ == "__main__":
    pytest.main()
