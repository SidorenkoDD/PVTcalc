"""Единый read-only renderer таблиц результатов."""

from collections.abc import Sequence

import dearpygui.dearpygui as dpg


def render_readonly_table(
    parent: int | str,
    columns: Sequence[object],
    rows: Sequence[Sequence[object]],
    *,
    scroll_x: bool = True,
    freeze_columns: int = 1,
) -> None:
    """Рисует результат как текстовую таблицу без редактируемых ячеек."""
    labels = [str(column) for column in columns]
    with dpg.table(
        parent=parent,
        header_row=True,
        borders_innerH=True,
        borders_outerH=True,
        borders_innerV=True,
        borders_outerV=True,
        resizable=True,
        scrollY=True,
        scrollX=scroll_x,
        height=-1,
        freeze_rows=1,
        freeze_columns=min(max(freeze_columns, 0), len(labels)),
    ):
        for label in labels:
            dpg.add_table_column(label=label)
        for row in rows:
            with dpg.table_row():
                values = list(row)
                for index in range(len(labels)):
                    value = values[index] if index < len(values) else ""
                    dpg.add_text(str(value))
