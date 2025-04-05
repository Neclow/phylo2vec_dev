"""Test conversion from matrix form to newick and back"""

import numpy as np
import pytest

from phylo2vec.tests.config import MIN_N_LEAVES, MAX_N_LEAVES, N_REPEATS
from phylo2vec.matrix import to_matrix, to_matrix_no_parents, to_newick
from phylo2vec.utils import remove_parent_labels, sample_matrix


@pytest.mark.parametrize("n_leaves", range(MIN_N_LEAVES, MAX_N_LEAVES + 1))
def test_m2newick2m(n_leaves):
    """Test that m to newick to converted_m via to_matrix leads to m == converted_m

    Parameters
    ----------
    n_leaves : int
        Number of leaves
    """
    for _ in range(N_REPEATS):
        m = sample_matrix(n_leaves)
        newick = to_newick(m)
        m2 = to_matrix(newick)
        assert np.allclose(m, m2, atol=1e-6)

        newick_no_parents = remove_parent_labels(newick)
        m3 = to_matrix_no_parents(newick_no_parents)
        assert np.allclose(m, m3, atol=1e-6)
