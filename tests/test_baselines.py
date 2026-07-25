import unittest

import numpy as np

from cisgrammar.baselines import PWMPresenceBaseline
from cisgrammar.data import Condition, GrammarRule, generate_matched_dataset
from cisgrammar.motifs import NANOG, POU5F1


class TestPWMPresenceBaseline(unittest.TestCase):
    def test_probabilities_are_finite(self):
        dataset = generate_matched_dataset(
            n_pairs=12,
            sequence_length=144,
            motif_a=POU5F1,
            motif_b=NANOG,
            rule=GrammarRule(),
            condition=Condition(gap_min=0, gap_max=50),
            anchor_min=6,
            anchor_max=12,
            seed=7,
            condition_name="baseline",
        )
        model = PWMPresenceBaseline(POU5F1, NANOG).fit(dataset.x)
        probabilities = model.predict_proba(dataset.x)
        self.assertEqual(probabilities.shape, (24,))
        self.assertTrue(np.all(np.isfinite(probabilities)))
        self.assertTrue(np.all((probabilities >= 0) & (probabilities <= 1)))


if __name__ == "__main__":
    unittest.main()
