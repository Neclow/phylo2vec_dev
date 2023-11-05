"""Base class for all optimisation methods in Phylo2Vec."""
import random

import numba as nb

from phylo2vec.datasets import load_fasta
from phylo2vec.utils import sample, seed_everything

MAX_SEED = 42


class BaseOptimizer:
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
        records = list(load_fasta(fasta_path))

        taxa_dict = self._make_taxa_dict(records)

        n_leaves = len(taxa_dict)

        v_init = sample(n_leaves)

        v_opt, losses = self._optimise(fasta_path, v_init, taxa_dict)

        return v_opt, losses

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
