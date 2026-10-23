"""Preservation and fail-closed checks for the notebook editing helper."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('notebook_cells', Path(__file__).parents[1] / 'tools/notebook_cells.py')
cells = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cells)


class NotebookCellsTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.path = self.root / 'example.ipynb'
        self.workspace = self.root / 'editing'
        self.notebook = {
            'nbformat': 4, 'nbformat_minor': 5, 'metadata': {'custom': 'keep'},
            'cells': [
                {'id': 'intro', 'cell_type': 'markdown', 'metadata': {}, 'source': ['Old title\n'],
                 'attachments': {'image.png': {'image/png': 'YWJj'}}},
                {'id': 'code', 'cell_type': 'code', 'metadata': {'tags': ['keep']},
                 'source': ['print("hello")\n'], 'execution_count': 7,
                 'outputs': [{'output_type': 'stream', 'name': 'stdout', 'text': ['hello\n']}]},
                {'id': 'raw', 'cell_type': 'raw', 'metadata': {}, 'source': 'raw text'},
            ],
        }
        self.original = json.dumps(self.notebook, indent=3).encode()
        self.path.write_bytes(self.original)
        cells.split(self.path, self.workspace)

    def test_noop_preserves_exact_bytes(self):
        cells.apply(self.workspace)
        self.assertEqual(self.path.read_bytes(), self.original)

    def test_markdown_preserves_every_other_field(self):
        (self.workspace / 'cells/0000.md').write_text('New title: θ\n', encoding='utf-8', newline='\n')
        cells.apply(self.workspace)
        expected = self.notebook
        expected['cells'][0]['source'] = ['New title: θ\n']
        self.assertEqual(json.loads(self.path.read_bytes()), expected)

    def test_code_requires_opt_in_and_clears_only_edited_outputs(self):
        (self.workspace / 'cells/0001.py').write_text('print("updated")\n', encoding='utf-8', newline='\n')
        with self.assertRaises(ValueError):
            cells.apply(self.workspace)
        self.assertEqual(self.path.read_bytes(), self.original)
        cells.apply(self.workspace, allow_code=True)
        expected = self.notebook
        expected['cells'][1].update(source=['print("updated")\n'], outputs=[], execution_count=None)
        self.assertEqual(json.loads(self.path.read_bytes()), expected)

    def test_missing_cell_refused(self):
        (self.workspace / 'cells/0001.py').unlink()
        with self.assertRaises(ValueError):
            cells.apply(self.workspace)
        self.assertEqual(self.path.read_bytes(), self.original)

    def test_extra_cell_refused(self):
        (self.workspace / 'cells/extra.py').write_text('pass')
        with self.assertRaises(ValueError):
            cells.apply(self.workspace)
        self.assertEqual(self.path.read_bytes(), self.original)

    def test_stale_original_refused(self):
        self.path.write_bytes(self.original + b'\n')
        with self.assertRaises(ValueError):
            cells.apply(self.workspace)
        self.assertEqual(self.path.read_bytes(), self.original + b'\n')

    def test_damaged_backup_refused(self):
        (self.workspace / 'original.ipynb').write_text('{}')
        with self.assertRaises(ValueError):
            cells.apply(self.workspace)
        self.assertEqual(self.path.read_bytes(), self.original)

    def test_existing_workspace_refused(self):
        with self.assertRaises(FileExistsError):
            cells.split(self.path, self.workspace)
        self.assertEqual((self.workspace / 'original.ipynb').read_bytes(), self.original)


if __name__ == '__main__':
    unittest.main()
