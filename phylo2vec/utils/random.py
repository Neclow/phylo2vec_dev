"""Random utility functions: sampling and seeding."""

import os
import random

import numba as nb
import numpy as np


@nb.njit
def sample(n_leaves, ordered=False):
    """Sample a random tree topology via Phylo2Vec

    Parameters
    ----------
    n_leaves : int
        Number of leaves
    ordered : bool, optional
        If True, sample an ordered tree, by default False

        True:
        v_i in {0, 1, ..., i} for i in (0, n_leaves-1)

        False:
        v_i in {0, 1, ..., 2*i} for i in (0, n_leaves-1)

    Returns
    -------
    numpy.ndarray
        Phylo2Vec vector
    """

    if ordered:
        v_list = [np.random.randint(0, i + 1) for i in range(n_leaves - 1)]
    else:
        v_list = [np.random.randint(0, 2 * i + 1) for i in range(n_leaves - 1)]
    return np.array(v_list, dtype=np.uint32)


BRANCH_LENGTH_DISTRIBUTIONS = {
    "exponential": np.random.exponential,
    "gamma": np.random.gamma,
    "uniform": np.random.uniform,
}


def sample_matrix(
    n_leaves, ordered=False, branch_length_distribution="uniform", **kwargs
):
    """Sample a random tree with branch lengths via Phylo2Vec

    By default, branch lengths are sampled uniformly in (0, 1)

    1st column: v (topology)
    2nd and 3rd columns: branch lengths of cherry $i$ in the ancestry matrix

    Parameters
    ----------
    n_leaves : int
        Number of leaves
    ordered : bool, optional
        If True, sample an ordered tree topology, by default False

        True:
        v_i in {0, 1, ..., i} for i in (0, n_leaves-1)

        False:
        v_i in {0, 1, ..., 2*i} for i in (0, n_leaves-1)
    branch_length_distribution : str, optional
        Distribution upon which branch lengths are distributed, by default 'uniform'
    kwargs :
        All optional arguments are passed to the branch length distribution

    Returns
    -------
    numpy.ndarray
        Phylo2Mat matrix
    """
    v = sample(n_leaves, ordered=ordered)

    bl_size = (v.shape[0], 2)

    if branch_length_distribution in BRANCH_LENGTH_DISTRIBUTIONS:
        dist = BRANCH_LENGTH_DISTRIBUTIONS[branch_length_distribution]
        bls = dist(size=bl_size, **kwargs)
    else:
        raise ValueError(
            f"`branch_length_distribution` must in be {BRANCH_LENGTH_DISTRIBUTIONS.keys()}. "
            f"Got {branch_length_distribution} instead."
        )

    m = np.concatenate([v[:, None], bls.astype(np.float16)], axis=1)

    return m


def seed_everything(seed):
    """Seed random, the Python hash seed, numpy

    Parameters
    ----------
    seed : int
        Random seed
    """
    random.seed(seed)
    os.environ["PYTHONHASHSEED"] = str(seed)
    np.random.seed(seed)
