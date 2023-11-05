"""Methods for hill-climbing optimisation."""
import numpy as np

from joblib import delayed, effective_n_jobs, Parallel

from phylo2vec.opt._base import BaseOptimizer
from phylo2vec.opt.hc._hc_losses import raxml_loss
from phylo2vec.utils.vector import reorder_v, reroot_at_random


class HillClimbingOptimizer(BaseOptimizer):
    def __init__(
        self,
        raxml_path="",
        tree_folder_path="",
        substitution_model="GTR",
        reorder_method="birth_death",
        random_seed=None,
        tol=0.001,
        patience=3,
        rounds=1,
        n_jobs=None,
        verbose=False,
    ):
        super().__init__(random_seed=random_seed)

        self.raxml_path = raxml_path
        self.tree_folder_path = tree_folder_path
        self.substitution_model = substitution_model
        self.reorder_method = reorder_method
        self.tol = tol
        self.rounds = rounds
        self.patience = patience
        self.n_jobs = effective_n_jobs(n_jobs)
        self.verbose = verbose

    def _optimise(self, fasta_path, v, taxa_dict):
        current_loss = raxml_loss(
            v=v,
            taxa_dict=taxa_dict,
            raxml_path=self.raxml_path,
            fasta_path=fasta_path,
            tree_folder_path=self.tree_folder_path,
            substitution_model=self.substitution_model,
        )

        wait = 0

        losses = [current_loss]

        while wait < self.patience:
            if self.verbose:
                print("Changing equivalences...")
            v_proposal = reroot_at_random(v)

            v_proposal, proposal_loss, taxa_dict = self._optimise_single(
                fasta_path, v.copy(), taxa_dict
            )

            v = v_proposal.copy()

            if proposal_loss - current_loss < self.tol:
                # Found a better loss so reset the patience counter
                current_loss = proposal_loss

                wait = 0
            else:
                # No drastic increase so increment the patience counter
                wait += 1

                if self.verbose:
                    print(f"No significantly better loss found {wait}/{self.patience}.")

            losses.append(current_loss)

        return v, losses

    def _optimise_single(self, fasta_path, v, taxa_dict):
        # Reorder v
        v_shuffled, taxa_dict = reorder_v(self.reorder_method, v, taxa_dict)

        # Get current loss
        current_best_loss = raxml_loss(
            v=v_shuffled,
            taxa_dict=taxa_dict,
            raxml_path=self.raxml_path,
            fasta_path=fasta_path,
            tree_folder_path=self.tree_folder_path,
            substitution_model=self.substitution_model,
            outfile="output.tree",
        )

        if self.verbose:
            print(f"Start optimise_single: {current_best_loss:.3f}")

        for _ in range(self.rounds):
            for i in reversed(range(1, len(v_shuffled))):
                # Calculate gradient for changes in row i
                # "gradient" here simply refers to a numerical gradient
                # between loss(v_current) and loss(v_proposal)
                proposal_grads, proposal_losses = self.grad_single(
                    fasta_path=fasta_path,
                    v_proposal=v_shuffled,
                    current_loss=current_best_loss,
                    taxa_dict=taxa_dict,
                    i=i,
                )

                # find index of max gradient
                grad_choice = proposal_grads.argmax(0)

                # Is there a positive gradient?
                if proposal_grads[grad_choice] > 0:
                    # Discrete gradient step
                    v_shuffled[i] = (
                        grad_choice + 1 if grad_choice >= v_shuffled[i] else grad_choice
                    )

                    # Reorder v
                    # v_shuffled, taxa_dict = reorder_v(v_shuffled, taxa_dict)

                    if self.verbose:
                        grad_propose = proposal_losses[grad_choice] - current_best_loss
                        print(
                            f"Loss: {proposal_losses[grad_choice]:.3f} (diff: {grad_propose:.3f})"
                        )

                    # Update best loss
                    current_best_loss = proposal_losses[grad_choice]

        if self.verbose:
            print(f"End optimise_single: {current_best_loss:.3f}")

        return v_shuffled, current_best_loss, taxa_dict

    def grad_single(self, fasta_path, v_proposal, current_loss, taxa_dict, i):
        v_copy = v_proposal.copy()

        def run(v_other, i, j):
            v_other[i] = j
            return raxml_loss(
                v=v_other,
                taxa_dict=taxa_dict,
                raxml_path=self.raxml_path,
                fasta_path=fasta_path,
                tree_folder_path=self.tree_folder_path,
                substitution_model=self.substitution_model,
                outfile=f"tree{i}{j}.tree",
            )

        proposal_losses = np.array(
            Parallel(n_jobs=self.n_jobs)(
                delayed(run)(v_copy, i, j)
                for j in range(2 * i + 1)
                if j != v_proposal[i]
            )
        )

        return current_loss - proposal_losses, proposal_losses
