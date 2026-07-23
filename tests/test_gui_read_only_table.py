"""Проверки выделения и TSV-копирования read-only таблиц."""

import dearpygui.dearpygui as dpg

from gui.view import read_only_table as table_view


def test_selection_tsv_contains_headers_and_rectangular_range():
    dpg.create_context()
    try:
        with dpg.window() as window:
            table_id = table_view.render_readonly_table(
                window, ["P", "Bo", "Rs"], [[300, 1.2, 100], [200, 1.1, 80]],
            )
        selection = table_view._TABLES[table_id]
        selection.anchor = (0, 1)
        selection.end = (1, 2)

        assert table_view._selection_tsv(table_id) == "Bo\tRs\n1.2\t100\n1.1\t80"
    finally:
        dpg.destroy_context()
