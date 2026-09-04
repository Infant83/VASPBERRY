"""Repository-level checks for files that must never be published."""

from pathlib import Path
import subprocess
import unittest


ROOT = Path(__file__).resolve().parents[1]


class RepositoryHygieneTests(unittest.TestCase):
    def test_no_potcar_file_is_tracked(self):
        output = subprocess.check_output(
            ["git", "ls-files", "-z"], cwd=ROOT
        ).decode("utf-8")
        tracked_potcars = sorted(
            path
            for path in output.split("\0")
            if path and Path(path).name.upper().startswith("POTCAR")
        )
        self.assertEqual(
            [],
            tracked_potcars,
            "VASP POTCAR data must not be tracked; remove: "
            + ", ".join(tracked_potcars),
        )


if __name__ == "__main__":
    unittest.main()
