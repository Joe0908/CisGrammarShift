import unittest

try:
    import torch
except ImportError:  # pragma: no cover - permits generator-only checks in minimal environments
    torch = None


@unittest.skipIf(torch is None, "PyTorch is not installed")
class TestModels(unittest.TestCase):
    def test_output_shapes(self):
        from cisgrammar.models import build_model

        x = torch.zeros(3, 128, 4)
        configurations = {
            "local_cnn": {"channels": 8, "kernel_size": 19},
            "dilated_cnn": {"channels": 8, "kernel_size": 19, "dilations": [1, 2, 4]},
            "transformer": {
                "embedding_dim": 16,
                "attention_heads": 4,
                "layers": 1,
                "feedforward_dim": 32,
                "dropout": 0.0,
            },
        }
        for name, configuration in configurations.items():
            with self.subTest(name=name):
                model = build_model(name, sequence_length=128, model_config=configuration)
                self.assertEqual(tuple(model(x).shape), (3,))


if __name__ == "__main__":
    unittest.main()
