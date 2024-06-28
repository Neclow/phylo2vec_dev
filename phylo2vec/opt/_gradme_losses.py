import jax.scipy as jscipy

from jax import jit


@jit
def get_edges_exp_log(W, rooted):
    raise NotImplementedError


@jit
def bme_loss_log(W, D, rooted):
    """Log version of the BME loss function

    Parameters
    ----------
    W : jax.numpy.array
        Ordered tree probability matrix
    D : jax.numpy.array
        Distance matrix
    rooted : bool
        True is the tree is rooted, otherwise False

    Returns
    -------
    float
        BME loss
    """
    E = get_edges_exp_log(W, rooted)

    return jscipy.special.logsumexp(E, b=D)
