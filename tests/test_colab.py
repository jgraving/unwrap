"""Unit tests for the Colab bootstrap helper (unwrap.colab)."""
import io
import sys
import tempfile
import types
import unittest
from importlib.metadata import PackageNotFoundError
from pathlib import Path
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


class ScriptsDirTests(unittest.TestCase):
    def test_targets_the_notebooks_dir_next_to_the_package(self):
        target = colab.scripts_dir()
        self.assertEqual(target.name, "notebooks")
        self.assertEqual(
            target.parent, Path(colab.__file__).resolve().parent.parent
        )


class ScriptFetchTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.target = Path(self.temp.name)
        self.opened = []
        # `import google.colab` resolves the parent package first, so
        # both levels of the namespace need stubs.
        self.parent = types.ModuleType("google")
        self.fake = types.ModuleType("google.colab")
        self.parent.colab = self.fake
        sys.modules["google"] = self.parent
        sys.modules["google.colab"] = self.fake

    def tearDown(self):
        del sys.modules["google.colab"]
        del sys.modules["google"]

    def fake_urlopen(self, url):
        self.opened.append(url)
        return io.BytesIO(b"# fetched script\n")

    def test_fetches_each_script_from_ref_into_scripts_dir(self):
        with mock.patch.object(colab, "scripts_dir", return_value=self.target), \
                mock.patch.object(colab.urllib.request, "urlopen",
                                  side_effect=self.fake_urlopen), \
                mock.patch("unwrap.colab.subprocess.check_call"):
            colab.install_colab(scripts=["a.py", "b.py"], ref="some-ref")
        self.assertEqual(self.opened, [
            "https://raw.githubusercontent.com/jgraving/unwrap/"
            "some-ref/notebooks/a.py",
            "https://raw.githubusercontent.com/jgraving/unwrap/"
            "some-ref/notebooks/b.py",
        ])
        self.assertEqual((self.target / "a.py").read_bytes(), b"# fetched script\n")
        self.assertEqual((self.target / "b.py").read_bytes(), b"# fetched script\n")


if __name__ == "__main__":
    unittest.main()
