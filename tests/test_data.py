import unittest

import numpy as np

from cisgrammar.data import Condition, GrammarRule, generate_matched_dataset, one_hot_encode
from cisgrammar.motifs import NANOG, POU5F1


class TestMatchedGenerator(unittest.TestCase):
    def setUp(self):
        self.rule = GrammarRule(period=10, allowed_phases=(0, 1, 9))
        self.dataset = generate_matched_dataset(
            n_pairs=24,
            sequence_length=160,
            motif_a=POU5F1,
            motif_b=NANOG,
            rule=self.rule,
            condition=Condition(
                gap_min=0,
                gap_max=60,
                background_gc=0.5,
                background_persistence=0.1,
                allow_reverse=True,
                motif_temperature=1.0,
            ),
            anchor_min=8,
            anchor_max=20,
            seed=123,
            condition_name="test",
        )

    def test_shapes_and_balance(self):
        self.assertEqual(self.dataset.x.shape, (48, 160, 4))
        self.assertEqual(self.dataset.motif_mask.shape, (48, 160))
        np.testing.assert_array_equal(np.bincount(self.dataset.y.astype(int)), [24, 24])
        np.testing.assert_allclose(self.dataset.x.sum(axis=-1), 1.0)

    def test_pair_invariants(self):
        for pair_id in np.unique(self.dataset.pair_ids):
            records = [record for record in self.dataset.records if record["pair_id"] == pair_id]
            self.assertEqual(sorted(record["label"] for record in records), [0, 1])
            self.assertEqual(records[0]["motif_a_instance"], records[1]["motif_a_instance"])
            self.assertEqual(records[0]["motif_b_instance"], records[1]["motif_b_instance"])
            self.assertEqual(records[0]["background_sha256"], records[1]["background_sha256"])
            for record in records:
                self.assertEqual(bool(record["label"]), self.rule.is_positive(record["gap"]))

    def test_reproducible_seed(self):
        repeated = generate_matched_dataset(
            n_pairs=24,
            sequence_length=160,
            motif_a=POU5F1,
            motif_b=NANOG,
            rule=self.rule,
            condition=Condition(
                gap_min=0,
                gap_max=60,
                background_gc=0.5,
                background_persistence=0.1,
                allow_reverse=True,
                motif_temperature=1.0,
            ),
            anchor_min=8,
            anchor_max=20,
            seed=123,
            condition_name="test",
        )
        np.testing.assert_array_equal(self.dataset.y, repeated.y)
        np.testing.assert_array_equal(self.dataset.pair_ids, repeated.pair_ids)

    def test_invalid_alphabet_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "unsupported base"):
            one_hot_encode(["ACNT"])


if __name__ == "__main__":
    unittest.main()
