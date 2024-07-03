import numba as nb
import numpy as np

from phylo2vec.base.to_newick import _get_ancestry
from phylo2vec.opt._base import BaseOptimizer


# class EdgeManager:
#     def __init__(self, v):
#         self.edges_r = self.get_rooted_edges(v)
#         self.edges_u = self.as_unrooted()

#     @nb.njit
#     def get_rooted_edges(self, v):
#         anc = _get_ancestry(v)

#         edges = np.zeros((2 * len(v), 2), dtype=np.int16)

#         for i in range(len(v)):
#             edges[2 * i] = anc[i, [0, 2]]
#             edges[2 * i + 1] = anc[i, [1, 2]]

#         return edges

#     @nb.njit
#     def as_unrooted(self):
#         edges_u = self.edges_r[:-1].copy()
#         edges_u[-2, :] = [self.edges_r.max() - 1, self.edges_r[-2:, :].min()]
#         return edges_u


class DiscreteOptimizer(BaseOptimizer):

    def _optimise(self, fasta_path, v, label_mapping):
        # calculate distance matrix
        raise NotImplementedError

    def grad(self, v, D):
        n_leaves = len(v) + 1
        edges_r = self.get_rooted_edges(v)
        edges_u = self.get_unrooted_edges(edges_r)
        edge_weights = self.get_edge_weights(n_leaves, D, edges_u)

        edge_locations_r = self.get_edge_locations(edges_r)
        edge_locations_u = self.get_edge_locations(edges_u)

    @staticmethod
    @nb.njit
    def get_rooted_edges(v: np.ndarray) -> np.ndarray:
        anc = _get_ancestry(v)

        edges = np.zeros((2 * len(v), 2), dtype=np.int16)

        for i in range(len(v)):
            edges[2 * i] = anc[i, [0, 2]]
            edges[2 * i + 1] = anc[i, [1, 2]]

        return edges

    @staticmethod
    @nb.njit
    def get_unrooted_edges(edges_r: np.ndarray) -> np.ndarray:
        edges_u = edges_r[:-1].copy()
        edges_u[-2, :] = [edges_r.max() - 1, edges_r[-2:, :].min()]
        return edges_u

    @staticmethod
    @nb.njit
    def get_edge_weights(n_leaves: int, D: np.ndarray, edges_u: list) -> np.ndarray:
        n_edges = 2 * n_leaves - 3

        edge_weights = np.zeros((n_edges, 2, n_leaves))

        edges_u_first = edges_u[:, 0]
        edges_u_last = edges_u[:, 1]

        # We begin by considering the subtrees not containing the root.
        for n in range(n_edges):
            child = edges_u[n, 1]

            # child is a leaf node
            if child < n_leaves:
                # This means that the subtree not containing the root is simply the child node.
                # Hence, we can set the subtree distances to simply be the corresponding row in D
                edge_weights[n, 0] = D[child]

            # child is an internal node
            else:
                # We must find the two subtrees connecting the child of the child node of edge i.
                # We must then average their distances.
                edges_grandchildren = edge_weights[edges_u_first == child, 0]

                edge_weights[n, 0] = edges_grandchildren.sum(0) / 2.0

        # We now consider the subtrees containing the root
        for n in range(n_edges - 2, -1, -1):
            # We process the edges in reverse order, so that we move further from the root (and
            # can hence calculate the edges iteratively as before).

            parent, child = edges_u[n]

            # Find edges connected to the parent node

            # This first set are edges that are further to the root as parent is their parent edge
            parent_filter = (edges_u_first == parent) & (edges_u_last != child)

            # This second set are edges that are further from the root
            child_filter = (edges_u_last == parent) & (edges_u_first != child)

            # Use previously calculated weights.

            # For these edges, we want the subtree not containing the root
            pweights = edge_weights[parent_filter, 0]

            # For these edges, we want the subtree containing the root
            cweights = edge_weights[child_filter, 1]

            # Add the total balanced contribution to edge weights
            # (note there have definitely been only two contributing edges)
            edge_weights[n, 1] = pweights.sum(0) + cweights.sum(0) / 2.0

        return edge_weights

    @staticmethod
    @nb.njit
    def get_edge_locations(edges):
        n_max = edges.max()
        edge_locations = np.zeros((n_max + 1, n_max + 1), dtype=edges.dtype)

        for e in edges:
            edge_locations[e[0], e[1]] = 1
            edge_locations[e[1], e[0]] = 1

        return edge_locations

    @staticmethod
    @nb.njit
    def make_subtree(n_leaves, edges_u, i):
        e = edges_u[i]

        tree_edges = [(e[0], e[1])]

        to_visit = e[0]

        visited = np.zeros((edges_u.max() + 1,), dtype=np.int16)

        visited[to_visit] = 1

        while True:
            next_visit = -1
            for j, new_e in enumerate(edges_u):
                if j == i:
                    continue
                if new_e[0] == to_visit and visited[new_e[1]] == 0:
                    next_visit = new_e[1]
                    tree_edges.insert(0, (new_e[1], new_e[0]))
                elif new_e[1] == to_visit and visited[new_e[0]] == 0:
                    next_visit = new_e[0]
                    tree_edges.insert(0, (new_e[0], new_e[1]))

            if next_visit > -1 and visited[next_visit] == 0:
                visited[next_visit] = 1
                to_visit = next_visit
            else:
                break

        tree_edges = np.array(tree_edges)

        tree_size = np.sum(tree_edges < n_leaves) - (tree_edges[-1, -1] < n_leaves)

        return tree_edges, tree_size + 1
