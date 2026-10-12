"""Sanity checks for the null-normalized statistics and coarse-graining.

Run from the project folder:  python -m unittest discover tests
"""
import unittest

import networkx as nx
import numpy as np

from sctransfer import normalized as N
from sctransfer.run_coarse_grain import coarse_grain


class NormalizedSignatureTest(unittest.TestCase):
    def test_random_graph_is_near_one(self):
        """An ER graph should be indistinguishable from its own nulls."""
        A = nx.to_numpy_array(nx.gnm_random_graph(150, 1200, seed=1))
        sig = N.normalized_signature(A, n_nulls=5, seed=1)
        for key in ["clustering_norm", "degree_cv_norm"]:
            self.assertAlmostEqual(sig[key], 1.0, delta=0.25, msg=key)

    def test_modular_graph_beats_random(self):
        A = nx.to_numpy_array(nx.planted_partition_graph(6, 25, 0.4, 0.01, seed=2))
        sig = N.normalized_signature(A, n_nulls=5, seed=2)
        self.assertGreater(sig["modularity_norm"], 1.5)
        self.assertGreater(sig["clustering_norm"], 1.5)

    def test_density_matching(self):
        rng = np.random.default_rng(3)
        W = rng.random((80, 80))
        A = N.binarize(W, density=0.1)
        self.assertAlmostEqual(N.density_of(A), 0.1, delta=0.01)

    def test_density_matching_with_tied_integer_weights(self):
        rng = np.random.default_rng(5)
        W = rng.integers(0, 3, size=(100, 100)).astype(float)  # many ties at 1 and 2
        A = N.binarize(W, density=0.05)
        self.assertAlmostEqual(N.density_of(A), 0.05, delta=0.002)

    def test_dense_graph_does_not_raise(self):
        A = nx.to_numpy_array(nx.gnp_random_graph(40, 0.9, seed=4))
        sig = N.normalized_signature(A, n_nulls=2, seed=4)
        self.assertTrue(np.isfinite(sig["modularity_norm"]))


class CoarseGrainTest(unittest.TestCase):
    def test_sums_between_group_weights(self):
        W = np.array([[0, 1, 2, 0],
                      [0, 0, 0, 3],
                      [0, 0, 0, 4],
                      [5, 0, 0, 0]], dtype=float)
        C, labels = coarse_grain(W, ["a", "a", "b", "b"])
        self.assertEqual(list(labels), ["a", "b"])
        np.testing.assert_array_equal(C, [[0, 2 + 3], [5, 0]])


if __name__ == "__main__":
    unittest.main()
