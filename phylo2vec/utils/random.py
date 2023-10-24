import random
import numpy as np


def sample(n_leaves):
    """Sample a random tree via Phylo2Vec

    Returns
    -------
    numpy.ndarray
        Phylo2Vec vector where v_i in {0, 1, ..., 2*i}
    """
    return np.array([random.randint(0, 2 * i) for i in range(n_leaves - 1)])
