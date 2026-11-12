"""Independent finite-size and data-integrity checks for the interactive layer."""
from fractions import Fraction
import hashlib
import unittest

import numpy as np

from explore import ROOT, experiment_metadata, load_manuscript, moment_table, run_experiment


class NumericalChecks(unittest.TestCase):
    def test_three_element_distribution(self):
        # At phi = 1/2, P(C_3 = 2) = 2/7 and P(C_3 = 3) = 5/7.
        result = run_experiment("0.5", 3, 10000, 73)
        self.assertAlmostEqual(float(result["mean"][3]), float(Fraction(19, 7)), places=12)
        for r, expected in ((2, Fraction(10, 49)), (3, Fraction(-30, 343)),
                            (4, Fraction(190, 2401))):
            self.assertAlmostEqual(float(result["moments"][r, 3]), float(expected), places=12)
        self.assertEqual(set(result["costs"]), {2, 3})
        count = np.count_nonzero(result["costs"] == 2)
        sd = np.sqrt(10000 * (2 / 7) * (5 / 7))
        self.assertLess(abs(count - 10000 * 2 / 7), 7 * sd)

    def test_original_numeric_files_are_unchanged(self):
        for line in (ROOT / "data/SHA256SUMS").read_text().splitlines():
            digest, relative = line.split("  ", 1)
            self.assertEqual(hashlib.sha256((ROOT / relative).read_bytes()).hexdigest(), digest,
                             msg=relative)

    def test_saved_statistics_and_metadata(self):
        for q in ("0.25", "0.5", "0.8"):
            result = load_manuscript(q)
            metadata = experiment_metadata(result)
            self.assertEqual(metadata["n"], 1000)
            self.assertEqual(metadata["repetitions"], 12000)
            self.assertEqual(metadata["max_centred_moment_order"], 6)
            self.assertTrue((result["costs"] <= 1000 * 999 // 2).all())
            table = moment_table(result)
            self.assertEqual(len(table), 1998)
            self.assertTrue(np.isfinite(table.to_numpy()).all())
            self.assertGreater(float(table["variance"].min()), 0)

    def test_seed_reproducibility_and_domain(self):
        first = run_experiment("0.8", 30, 500, 17)
        second = run_experiment("0.8", 30, 500, 17)
        np.testing.assert_array_equal(first["costs"], second["costs"])
        for phi in ("0", "1", "nan", "inf"):
            with self.assertRaises(ValueError):
                run_experiment(phi)
        with self.assertRaises(ValueError):
            run_experiment("0.5", n=2)


if __name__ == "__main__":
    unittest.main()
