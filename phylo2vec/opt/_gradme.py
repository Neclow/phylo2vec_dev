from phylo2vec.opt._base import BaseOptimizer


class GradMEOptimizer(BaseOptimizer):

    def _optimise(self, fasta_path, v, label_mapping):
        raise NotImplementedError
