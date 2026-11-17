import ast
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "XRD_Analyzer_GUI.py"


def load_normalizer():
    tree = ast.parse(SOURCE.read_text(encoding="utf-8"))
    wanted = []
    for node in tree.body:
        if isinstance(node, ast.Assign):
            names = {t.id for t in node.targets if isinstance(t, ast.Name)}
            if "DEFAULT_SEPARATOR_COLOR" in names:
                wanted.append(node)
        elif isinstance(node, ast.FunctionDef) and node.name == "normalizuj_separator":
            wanted.append(node)
    ns = {}
    exec(compile(ast.Module(body=wanted, type_ignores=[]), str(SOURCE), "exec"), ns, ns)
    return ns["normalizuj_separator"]


normalizuj_separator = load_normalizer()


class SeparatorMigrationTests(unittest.TestCase):
    def test_old_single_label_moves_above_line(self):
        sep = normalizuj_separator({"position": 4, "label": "Materiały odniesienia", "offset_pct": 12}, 8)
        self.assertEqual(sep["top_label"], "Materiały odniesienia")
        self.assertEqual(sep["bottom_label"], "")
        self.assertEqual(sep["offset_pct"], 12.0)
        self.assertEqual(sep["color"], "#808080")
        self.assertTrue(sep["labels_visible"])

    def test_new_schema_preserves_two_labels(self):
        src = {
            "position": 4,
            "top_label": "Materiały odniesienia",
            "bottom_label": "Materiały badane",
            "offset_pct": -18,
            "color": "#444444",
            "labels_visible": False,
        }
        self.assertEqual(normalizuj_separator(src, 8), src)

    def test_bounds_and_invalid_values(self):
        self.assertEqual(normalizuj_separator({"position": -7, "offset_pct": 500}, 8)["position"], 0)
        self.assertEqual(normalizuj_separator({"position": 99, "offset_pct": 500}, 8)["position"], 8)
        self.assertEqual(normalizuj_separator({"offset_pct": 500}, 8)["offset_pct"], 48.0)
        self.assertEqual(normalizuj_separator({"offset_pct": -500}, 8)["offset_pct"], -48.0)
        self.assertEqual(normalizuj_separator({"offset_pct": "bad"}, 8)["offset_pct"], 0.0)


if __name__ == "__main__":
    unittest.main()
