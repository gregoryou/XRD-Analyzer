import ast
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "XRD_Analyzer_GUI.py"
TEXT = SOURCE.read_text(encoding="utf-8")
TREE = ast.parse(TEXT)


def method_source(class_name, method_name):
    cls = next(n for n in TREE.body if isinstance(n, ast.ClassDef) and n.name == class_name)
    fn = next(n for n in cls.body if isinstance(n, ast.FunctionDef) and n.name == method_name)
    lines = TEXT.splitlines()
    return "\n".join(lines[fn.lineno-1:fn.end_lineno])


class SeparatorRenderingStaticTests(unittest.TestCase):
    def test_offset_points_and_zorder(self):
        src = method_source("PlotCreatorDialog", "_draw_separator")
        self.assertIn('textcoords="offset points"', src)
        self.assertIn('zorder=0.4', src)
        self.assertIn('zorder=7', src)
        self.assertIn('labels_visible', src)

    def test_collision_solver_checks_peak_labels_and_traces(self):
        src = method_source("PlotCreatorDialog", "_resolve_separator_collisions")
        self.assertIn('peak_boxes', src)
        self.assertIn('_bbox_koliduje_z_linia', src)
        self.assertIn('x_candidates', src)
        self.assertIn('occupied', src)

    def test_multiple_same_boundary_separators_get_render_only_shift(self):
        src = method_source("PlotCreatorDialog", "redraw")
        self.assertIn('_data_dy_from_points', src)
        self.assertIn('groups.setdefault', src)
        self.assertIn('_draw_separator', src)


if __name__ == "__main__":
    unittest.main()
