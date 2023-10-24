import numpy as np


def check_v(v):
    k = len(v)

    v_max = 2 * np.arange(k)

    assert np.all((0 <= v) & (v <= v_max)), print(v, v >= 0, v <= v_max)
