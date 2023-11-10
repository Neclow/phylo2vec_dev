"""Base class for all optimisation methods in Phylo2Vec."""
import random

import numba as nb

from phylo2vec.datasets import read_fasta
from phylo2vec.utils import sample, seed_everything

MAX_SEED = 42


class BaseOptimizer:
    """
    Base class for all phylo2vec-based optimizers

    Parameters
    ----------
    random_seed : int, optional
        Random seed, by default None
    """

    def __init__(self, random_seed=None):
        self.random_seed = (
            random.randint(0, MAX_SEED) if random_seed is None else random_seed
        )
        seed_everything(self.random_seed)

    @staticmethod
    def _make_taxa_dict(records):
        taxa_dict = nb.typed.Dict.empty(
            key_type=nb.types.int64, value_type=nb.types.unicode_type
        )

        for i, r in enumerate(records):
            taxa_dict[i] = r.id.replace(" ", ".")

        return taxa_dict

    def fit(self, fasta_path):
        """Fit an optimizer to a fasta file

        Parameters
        ----------
        fasta_path : str
            Path to fasta file

        Returns
        -------
        v_opt : numpy.ndarray
            Optimized phylo2vec vector
        taxa_dict : dict[int, str]
            Mappping of leaf labels (integer) to taxa
        losses : array-like
            List/Array of collected losses
        """
        # TODO: figure out how to change this when the user selects load fasta
        # Probably an boolean "preloaded"
        # If True, change the fasta path to include the module name
        records = list(read_fasta(fasta_path))

        taxa_dict = self._make_taxa_dict(records)

        n_leaves = len(taxa_dict)

        v_init = sample(n_leaves)

        v_opt, taxa_dict, losses = self._optimise(fasta_path, v_init, taxa_dict)

        return v_opt, taxa_dict, losses

    def _optimise(self, fasta_path, v, taxa_dict):
        raise NotImplementedError

    def __repr__(self):
        # TODO: maybe something like sklearn pprint
        # https://github.com/scikit-learn/scikit-learn/blob/093e0cf14aff026cca6097e8c42f83b735d26358/sklearn/utils/_pprint.py#L116
        format_string = f"{self.__class__.__name__}("

        for item in vars(self):
            format_string += "\n"
            # TODO: pprint if dict?
            format_string += f"\t{item}={repr(self.__getattribute__(item))},"

        format_string = format_string[:-1] + "\n)"

        return format_string
