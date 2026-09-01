"""Optional integration test for a private full-mesh 1H-MoS2 WAVECAR.

Set ``VASPBERRY_MOS2_WAVECAR`` to enable it.  The fixture itself and any hash,
absolute path, or derived output are intentionally absent from the repository.
"""

import importlib.util
import os
import sys
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "wavecar_z2_integration", ROOT / "tools" / "wavecar_z2.py"
)
Z2 = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
sys.modules[SPEC.name] = Z2
SPEC.loader.exec_module(Z2)


FIXTURE = os.environ.get("VASPBERRY_MOS2_WAVECAR")


@unittest.skipUnless(FIXTURE, "set VASPBERRY_MOS2_WAVECAR for the private integration fixture")
class MoS2Z2IntegrationTest(unittest.TestCase):
    def test_12x12_reader_and_validation_contract(self):
        """Exercise the private fixture without baking in an unverified parity.

        The fixture was not available when version 1.1.0 was finalized, so a
        candidate value or a particular failed convergence guard must not be
        asserted here until the private full-mesh WAVECAR is rerun with the
        released implementation.
        """

        wavecar = Z2.Wavecar(Path(FIXTURE), spinor_components=2)
        result = Z2.analyze_wavecar(
            wavecar,
            nx=12,
            ny=12,
            occupied=18,
            axes=(0, 1),
            thresholds=Z2.Thresholds(),
        )
        self.assertIn(result["candidate_z2"], (0, 1, None))
        if result["valid"]:
            self.assertEqual(result["z2"], result["candidate_z2"])
            self.assertFalse(result["failed_guards"])
        else:
            self.assertIsNone(result["z2"])
            self.assertTrue(result["failed_guards"])
        self.assertLess(result["metrics"]["max_projector_tr_residual"], 1.0e-6)
        self.assertLess(abs(result["metrics"]["chern"]), 1.0e-4)
        self.assertGreater(result["metrics"]["direct_gap_ev"], 1.0)
        self.assertAlmostEqual(
            result["axes"]["x"]["candidate_z2"],
            result["axes"]["y"]["candidate_z2"],
        )


if __name__ == "__main__":
    unittest.main()
