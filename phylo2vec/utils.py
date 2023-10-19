import random

import numpy as np


def check_v(v):
    k = len(v)

    v_max = 2 * np.arange(k)

    assert np.all((0 <= v) & (v <= v_max)), print(v, v >= 0, v <= v_max)


def sample(k):
    """Sample a random tree via Phylo2Vec

    Returns
    -------
    numpy.ndarray
        Phylo2Vec vector where v_i in {0, 1, ..., 2*i}
    """
    return np.array([random.randint(0, 2 * i) for i in range(k)])
