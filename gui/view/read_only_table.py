"""Read-only tables with spreadsheet-like selection and clipboard support."""

from collections.abc import Sequence
from dataclasses import dataclass, field

import dearpygui.dearpygui as dpg

from gui.view.table_clipboard import table_to_tsv


@dataclass
class _TableSelection:
    columns: list[str]
    rows: list[list[str]]
    cells: dict[tuple[int, int], int] = field(default_factory=dict)
    anchor: tuple[int, int] | None = None
    end: tuple[int, int] | None = None


_TABLES: dict[int, _TableSelection] = {}
_ACTIVE_TABLE: int | None = None
_KEY_REGISTRY: int | None = None
_CELL_THEME: int | None = None
_SELECTED_CELL_THEME: int | None = None


def _is_dpg_type(item_id: int | None, expected: str) -> bool:
    """DPG item ids can be reused after a destroyed test/application context."""
    if item_id is None or not dpg.does_item_exist(item_id):
        return False
    try:
        return dpg.get_item_info(item_id).get("type") == expected
    except Exception:  # noqa: BLE001 - stale ids are not queryable on all DPG builds
        return False


def _cell_theme(*, selected: bool) -> int:
    global _CELL_THEME, _SELECTED_CELL_THEME
    current = _SELECTED_CELL_THEME if selected else _CELL_THEME
    if current is not None and _is_dpg_type(current, "mvAppItemType::mvTheme"):
        return current
    with dpg.theme() as theme_id:
        with dpg.theme_component(dpg.mvButton):
            dpg.add_theme_color(dpg.mvThemeCol_Button,
                                (205, 219, 232, 255) if selected else (0, 0, 0, 0))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonHovered, (213, 225, 237, 255))
            dpg.add_theme_color(dpg.mvThemeCol_ButtonActive, (195, 212, 228, 255))
            dpg.add_theme_color(dpg.mvThemeCol_Border,
                                (112, 145, 177, 255) if selected else (0, 0, 0, 0))
            dpg.add_theme_style(dpg.mvStyleVar_FrameBorderSize, 1 if selected else 0)
            dpg.add_theme_style(dpg.mvStyleVar_ButtonTextAlign, 0.02, 0.5)
    if selected:
        _SELECTED_CELL_THEME = theme_id
    else:
        _CELL_THEME = theme_id
    return theme_id


def _selection_bounds(selection: _TableSelection) -> tuple[int, int, int, int] | None:
    if selection.anchor is None or selection.end is None:
        return None
    row_a, column_a = selection.anchor
    row_b, column_b = selection.end
    return min(row_a, row_b), max(row_a, row_b), min(column_a, column_b), max(column_a, column_b)


def _is_selected(selection: _TableSelection, row: int, column: int) -> bool:
    bounds = _selection_bounds(selection)
    if bounds is None:
        return False
    row_min, row_max, column_min, column_max = bounds
    return row_min <= row <= row_max and column_min <= column <= column_max


def _refresh_selection(table_id: int) -> None:
    selection = _TABLES.get(table_id)
    if selection is None:
        return
    for (row, column), item_id in selection.cells.items():
        if dpg.does_item_exist(item_id):
            dpg.bind_item_theme(item_id, _cell_theme(
                selected=_is_selected(selection, row, column)))


def _on_cell_click(sender, app_data, user_data) -> None:
    del sender, app_data
    global _ACTIVE_TABLE
    table_id, row, column = user_data
    selection = _TABLES.get(table_id)
    if selection is None:
        return
    _ACTIVE_TABLE = table_id
    if (dpg.is_key_down(dpg.mvKey_ModShift)
            or dpg.is_key_down(dpg.mvKey_LShift)
            or dpg.is_key_down(dpg.mvKey_RShift)):
        selection.end = (row, column)
    else:
        selection.anchor = (row, column)
        selection.end = (row, column)
    _refresh_selection(table_id)


def _selection_tsv(table_id: int) -> str | None:
    selection = _TABLES.get(table_id)
    bounds = _selection_bounds(selection) if selection else None
    if selection is None or bounds is None:
        return None
    row_min, row_max, column_min, column_max = bounds
    columns = selection.columns[column_min:column_max + 1]
    rows = [row[column_min:column_max + 1]
            for row in selection.rows[row_min:row_max + 1]]
    return table_to_tsv(columns, rows)


def _copy_selection(table_id: int) -> bool:
    text = _selection_tsv(table_id)
    if text is None:
        return False
    try:
        dpg.set_clipboard_text(text)
    except Exception:  # noqa: BLE001 - clipboard may be unavailable in a remote session
        return False
    return True


def _on_copy_selection(sender, app_data, user_data) -> None:
    del sender, app_data
    _copy_selection(int(user_data))


def _on_copy_table(sender, app_data, user_data) -> None:
    del sender, app_data
    selection = _TABLES.get(int(user_data))
    if selection is not None:
        dpg.set_clipboard_text(table_to_tsv(selection.columns, selection.rows))


def _on_key_copy(sender, app_data) -> None:
    del sender, app_data
    if not (dpg.is_key_down(dpg.mvKey_ModCtrl)
            or dpg.is_key_down(dpg.mvKey_LControl)
            or dpg.is_key_down(dpg.mvKey_RControl)):
        return
    if _ACTIVE_TABLE is not None and dpg.does_item_exist(_ACTIVE_TABLE):
        _copy_selection(_ACTIVE_TABLE)


def _ensure_copy_handler() -> None:
    global _KEY_REGISTRY
    if _is_dpg_type(_KEY_REGISTRY, "mvAppItemType::mvHandlerRegistry"):
        return
    with dpg.handler_registry() as registry:
        dpg.add_key_press_handler(dpg.mvKey_C, callback=_on_key_copy)
    _KEY_REGISTRY = registry


def render_readonly_table(
    parent: int | str,
    columns: Sequence[object],
    rows: Sequence[Sequence[object]],
    *,
    scroll_x: bool = True,
    freeze_columns: int = 1,
) -> int:
    """Renders a selectable table; click then Shift+click selects a rectangle."""
    _ensure_copy_handler()
    labels = [str(column) for column in columns]
    values = [
        [str(row[index]) if index < len(row) and row[index] is not None else ""
         for index in range(len(labels))]
        for row in rows
    ]
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
    ) as table_id:
        _TABLES[table_id] = _TableSelection(labels, values)
        for label in labels:
            dpg.add_table_column(label=label)
        for row_index, row in enumerate(values):
            with dpg.table_row():
                for column_index, value in enumerate(row):
                    cell_id = dpg.add_button(
                        label=value or " ", width=-1,
                        user_data=(table_id, row_index, column_index),
                        callback=_on_cell_click,
                    )
                    _TABLES[table_id].cells[(row_index, column_index)] = cell_id
                    dpg.bind_item_theme(cell_id, _cell_theme(selected=False))
                    with dpg.popup(cell_id, mousebutton=dpg.mvMouseButton_Right):
                        dpg.add_menu_item(label="Copy selection", user_data=table_id,
                                          callback=_on_copy_selection)
                        dpg.add_menu_item(label="Copy table", user_data=table_id,
                                          callback=_on_copy_table)
    return table_id
