import unittest

import numpy as np

from cisgrammar.metrics import binary_metrics, paired_metrics, select_threshold


class TestMetrics(unittest.TestCase):
    def test_perfect_predictions(self):
        y = np.array([0, 0, 1, 1])
        probabilities = np.array([0.1, 0.2, 0.8, 0.9])
        threshold = select_threshold(y, probabilities)
        metrics = binary_metrics(y, probabilities, threshold)
        self.assertEqual(metrics["auroc"], 1.0)
        self.assertEqual(metrics["balanced_accuracy"], 1.0)
        self.assertEqual(metrics["mcc"], 1.0)

    def test_paired_metrics(self):
        y = np.array([1, 0, 0, 1])
        probabilities = np.array([0.9, 0.2, 0.7, 0.4])
        pairs = np.array(["a", "a", "b", "b"])
        metrics = paired_metrics(y, probabilities, pairs)
        self.assertEqual(metrics["pairwise_accuracy"], 0.5)
        self.assertAlmostEqual(metrics["counterfactual_delta"], 0.2)

    def test_malformed_pair_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "one-positive"):
            paired_metrics(
                np.array([1, 1]),
                np.array([0.6, 0.7]),
                np.array(["a", "a"]),
            )


if __name__ == "__main__":
    unittest.main()
