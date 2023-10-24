import random
import re

import numpy as np


LEFT_NODE_PATTERN = re.compile(r"\(\b\d+\b")
RIGHT_NODE_PATTERN = re.compile(r",\b\d+\b")
ANNOTATION_PATTERN = re.compile(r":\d+(\.\d+)?")


def check_v(v):
    k = len(v)

    v_max = 2 * np.arange(k)

    assert np.all((0 <= v) & (v <= v_max)), print(v, v >= 0, v <= v_max)


def find_num_leaves(newick):
    return len(re.findall(LEFT_NODE_PATTERN, newick)) + len(
        re.findall(RIGHT_NODE_PATTERN, newick)
    )


def remove_annotations_from_newick(newick):
    return re.sub(ANNOTATION_PATTERN, "", newick)


def sample(n_leaves):
    """Sample a random tree via Phylo2Vec

    Returns
    -------
    numpy.ndarray
        Phylo2Vec vector where v_i in {0, 1, ..., 2*i}
    """
    return np.array([random.randint(0, 2 * i) for i in range(n_leaves - 1)])
