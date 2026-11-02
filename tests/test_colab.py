"""Unit tests for the Colab bootstrap helper (unwrap.colab)."""
import sys
import types
import unittest
from importlib.metadata import PackageNotFoundError
from unittest import mock

from unwrap import colab


class NotebookExtraRequirementsTests(unittest.TestCase):
    def test_returns_the_notebook_extra_without_markers(self):
        requirements = colab.notebook_extra_requirements()
        self.assertGreater(len(requirements), 0)
        self.assertTrue(any(r.startswith("seaborn") for r in requirements))
        for req in requirements:
            self.assertNotIn(";", req)
            self.assertNotIn("extra", req)

    def test_missing_distribution(self):
        with mock.patch.object(
            colab, "requires", side_effect=PackageNotFoundError("unwrap")
        ):
            self.assertEqual(colab.notebook_extra_requirements(), [])


class InstallColabTests(unittest.TestCase):
    def test_noop_outside_colab(self):
        with mock.patch("unwrap.colab.subprocess.check_call") as call:
            colab.install_colab()
        call.assert_not_called()

    def test_installs_extra_and_per_notebook_extras_in_colab(self):
        # `import google.colab` resolves the parent package first, so
        # both levels of the namespace need stubs.
        parent = types.ModuleType("google")
        fake = types.ModuleType("google.colab")
        parent.colab = fake
        sys.modules["google"] = parent
        sys.modules["google.colab"] = fake
        try:
            with mock.patch("unwrap.colab.subprocess.check_call") as call:
                colab.install_colab(extras=["osmnx"])
            argv = call.call_args[0][0]
            self.assertEqual(argv[:5], [sys.executable, "-m", "pip", "install", "-q"])
            self.assertTrue(any(r.startswith("seaborn") for r in argv))
            self.assertEqual(argv[-1], "osmnx")
        finally:
            del sys.modules["google.colab"]
            del sys.modules["google"]


if __name__ == "__main__":
    unittest.main()
