# These tests were generated with Claude (Anthropic) and then reviewed, edited and run by the maintainers.
"""
Tests that keep tileseq_demux a separate preprocessor: no module in TileSeqDemux imports TileSeqMut, and no module in
TileSeqMut imports TileSeqDemux. The two programs meet only at the parameter sheet and the output folder.

Run with:  python -m pytest tests/test_demux_decoupling.py -v     (from the repository root)
Without pytest: python -m unittest discover -s tests -p "test_demux_*.py"
"""
import ast
import os
import subprocess
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def python_files(package):
    folder = os.path.join(ROOT, package)
    for dirpath, _dirs, files in os.walk(folder):
        if "__pycache__" in dirpath:
            continue
        for name in files:
            if name.endswith(".py"):
                yield os.path.join(dirpath, name)


def imported_packages(path):
    """Top-level names imported by one file."""
    with open(path, encoding="utf-8") as fh:
        tree = ast.parse(fh.read(), filename=path)
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            names.update(alias.name.split(".")[0] for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            names.add(node.module.split(".")[0])
    return names


class BoundaryTests(unittest.TestCase):
    def assert_no_import(self, package, forbidden):
        files = list(python_files(package))
        self.assertTrue(files, "no Python files found under %s; is the test run from the repository root?" % package)
        offenders = [os.path.relpath(f, ROOT) for f in files if forbidden in imported_packages(f)]
        self.assertEqual(offenders, [], "%s must not import %s" % (package, forbidden))

    def test_demux_does_not_import_tileseq_mut(self):
        self.assert_no_import("TileSeqDemux", "TileSeqMut")

    def test_tileseq_mut_does_not_import_demux(self):
        self.assert_no_import("TileSeqMut", "TileSeqDemux")

    def test_importing_demux_does_not_load_tileseq_mut(self):
        """In a fresh interpreter, so other tests cannot have loaded TileSeqMut already."""
        code = ("import sys; import TileSeqDemux.sheet, TileSeqDemux.validation; "
                "sys.exit(1 if any(m == 'TileSeqMut' or m.startswith('TileSeqMut.') for m in sys.modules) else 0)")
        result = subprocess.run([sys.executable, "-c", code], cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        self.assertEqual(result.returncode, 0, result.stdout)

    def test_demux_package_has_its_own_init(self):
        self.assertTrue(os.path.isfile(os.path.join(ROOT, "TileSeqDemux", "__init__.py")))


if __name__ == "__main__":
    unittest.main()