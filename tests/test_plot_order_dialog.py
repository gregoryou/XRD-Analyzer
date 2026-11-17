import ast
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT = (ROOT / "XRD_Analyzer_GUI.py").read_text(encoding="utf-8")
TREE = ast.parse(TEXT)


class PlotOrderStaticTests(unittest.TestCase):
    def test_dialog_has_two_label_fields_and_visibility(self):
        self.assertIn("Nazwa próbki / Napis nad linią", TEXT)
        self.assertIn("Plik źródłowy / Napis pod linią", TEXT)
        self.assertIn("labels_visible", TEXT)
        self.assertIn("Wybierz kolor linii podziału", TEXT)

    def test_project_version_is_13(self):
        main = next(n for n in TREE.body if isinstance(n, ast.ClassDef) and n.name == "MainWindow")
        value = None
        for n in main.body:
            if isinstance(n, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "PROJECT_VERSION" for t in n.targets):
                value = ast.literal_eval(n.value)
        self.assertEqual(value, 13)

    def test_legacy_group_separator_migration_is_kept(self):
        self.assertIn('state.get("group_separator", False)', TEXT)
        self.assertIn('state.get("top_group_label", "Reference materials")', TEXT)
        self.assertIn('state.get("bottom_group_label", "Dissimilar welded joints")', TEXT)


if __name__ == "__main__":
    unittest.main()
