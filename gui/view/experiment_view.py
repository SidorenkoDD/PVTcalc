"""Вкладки PVT-экспериментов, таблицы и графики результатов."""

import math

import dearpygui.dearpygui as dpg

from gui.app_state import NodeStatus
from gui.services import clipboard_service
from gui.services import experiment_service as exp_svc
from gui.view.contracts import ContextBoundView


class ExperimentViewMixin(ContextBoundView):
    """Рендер и callbacks CCE/DLE/Separator."""

    _exp_input_ids: dict[str, dict[str, int]]
    _exp_chart_holder: dict[str, int]
    _lab_data_holder: dict[str, int]
    _lab_data_controls: dict[str, tuple[int, int, int]]
    _lab_focus_theme_id: int | str | None
    _lab_active_cell: tuple[str, int, int] | None
    _lab_cell_ids: dict[tuple[str, int, int], int]
    _lab_navigation_registry_id: int | None

    def _render_experiment_tab(self, parent, node) -> None:
        nid = node.node_id
        p = node.params
        kind = p.get("kind", "cce")
        meta = exp_svc.EXPERIMENT_TYPES.get(kind, {})
        running = node.status is NodeStatus.RUNNING

        dpg.add_text(f"{meta.get('title', kind.upper())}", parent=parent)
        ids: dict = {}
        ids["pressures"] = dpg.add_input_text(
            label="Pressure stages, bar (comma-separated)", width=460, parent=parent,
            default_value=", ".join(self._g(x) for x in p.get("pressures", [])))
        with dpg.group(horizontal=True, parent=parent):
            ids["T_c"] = dpg.add_input_float(label="T, C", width=140, step=0,
                                             default_value=float(p.get("T_c", 100.0)))
            if meta.get("needs_p_res"):
                ids["P_res"] = dpg.add_input_float(
                    label="P res, bar", width=140, step=0,
                    default_value=float(p.get("P_res", 400.0)))
        if meta.get("needs_stage_temps"):
            ids["stage_temps"] = dpg.add_input_text(
                label="Stage T, C (comma-separated, same count)", width=460,
                parent=parent,
                default_value=", ".join(self._g(x) for x in p.get("stage_temps_c", [])))
        self._exp_input_ids[nid] = ids

        with dpg.group(horizontal=True, parent=parent):
            if running:
                dpg.add_loading_indicator(style=1, radius=2.0)
                dpg.add_button(label="Cancel", user_data=nid,
                               callback=self._on_flash_cancel)
            else:
                dpg.add_button(label="Run experiment", user_data=nid,
                               callback=self._on_experiment_run)
        if running:
            dpg.add_text("Running experiment (this may take a few seconds)...",
                         parent=parent)
        if node.status is NodeStatus.STALE and node.error:
            err = dpg.add_text(f"Error: {node.error}", parent=parent)
            dpg.bind_item_theme(err, self._theme_stale())
        elif node.status is NodeStatus.STALE and node.result is not None:
            st = dpg.add_text("Result is stale (composition changed) - re-run.",
                              parent=parent)
            dpg.bind_item_theme(st, self._theme_stale())

        self._render_lab_data(node, parent)

        if node.result is not None:
            dpg.add_separator(parent=parent)
            with dpg.tab_bar(parent=parent):
                with dpg.tab(label="Main") as t:
                    self._render_exp_table(node.result, t,
                                           node.result.get("main_columns", []))
                with dpg.tab(label="All data") as t:
                    self._render_exp_table(node.result, t, node.result["columns"])
                with dpg.tab(label="Composition by stage") as t:
                    self._render_exp_composition(node.result, t)
                with dpg.tab(label="Chart") as t:
                    self._render_exp_chart(node.result, t, nid)

    def _lab_columns(self, node) -> list[str]:
        kind = node.params.get("kind", "cce")
        meta = exp_svc.EXPERIMENT_TYPES.get(kind, {})
        return list(meta.get("lab_columns", ["pressure"]))

    def _lab_rows(self, node, columns: list[str]) -> list[list[float | None]]:
        data = node.params.get("lab_data")
        raw_rows = data.get("rows") if isinstance(data, dict) else None
        if not isinstance(raw_rows, list):
            return []
        rows: list[list[float | None]] = []
        for raw in raw_rows:
            if not isinstance(raw, list):
                continue
            rows.append([
                value if isinstance(value, (int, float))
                and not isinstance(value, bool) and math.isfinite(float(value))
                else None
                for value in (raw[:len(columns)] + [None] * len(columns))[:len(columns)]
            ])
        return rows

    def _render_lab_data(self, node, parent) -> None:
        """Сворачиваемый ввод измерений, привязанный к одному Experiment-узлу."""
        columns = self._lab_columns(node)
        rows = self._lab_rows(node, columns)
        with dpg.collapsing_header(label="Lab Data (measured)",
                                  default_open=True, parent=parent) as header:
            dpg.add_text(
                "Enter measured points for this experiment. Blank cells are ignored.",
                parent=header,
            )
            with dpg.group(horizontal=True, parent=header):
                dpg.add_button(label="Add point", user_data=node.node_id,
                               callback=self._on_lab_add_row)
                dpg.add_button(label="Paste from Excel",
                               user_data=node.node_id,
                               callback=self._on_lab_paste,
                               enabled=False)
                remove_id = dpg.add_button(
                    label="Remove last", user_data=node.node_id,
                    callback=self._on_lab_remove_row, enabled=bool(rows))
                clear_id = dpg.add_button(
                    label="Clear", user_data=node.node_id,
                    callback=self._on_lab_clear, enabled=bool(rows))
                count_id = dpg.add_text(f"{len(rows)} point(s)", parent=header)
                self._lab_data_controls[node.node_id] = (
                    remove_id, clear_id, count_id)
            holder = dpg.add_group(parent=header)
            self._lab_data_holder[node.node_id] = holder
            dpg.add_text(
                "Excel paste: tab-separated columns; a header row is detected "
                "automatically. Paste starts at the focused cell (or the first "
                "cell); Ctrl+V and the button both support column paste; arrow "
                "keys move between cells.",
                parent=header,
            )
            self._render_lab_data_table(node, holder, columns, rows)

    def _render_lab_data_table(self, node, holder, columns, rows) -> None:
        self._ensure_lab_navigation_handlers()
        self._lab_cell_ids = {
            key: cell_id for key, cell_id in self._lab_cell_ids.items()
            if key[0] != node.node_id
        }
        if dpg.does_item_exist(holder):
            dpg.delete_item(holder, children_only=True)
        if not rows:
            dpg.add_text("No measured points yet.", parent=holder)
            return
        with dpg.table(parent=holder, header_row=True,
                       borders_innerH=True, borders_outerH=True,
                       borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True,
                       height=190, freeze_rows=1):
            dpg.add_table_column(label="#", width_fixed=True, width=35)
            for column in columns:
                dpg.add_table_column(label=column)
            for row_index, row in enumerate(rows):
                with dpg.table_row():
                    dpg.add_text(str(row_index + 1))
                    for column_index, value in enumerate(row):
                        cell_id = dpg.add_input_text(
                            default_value="" if value is None else self._g(value),
                            width=115,
                            user_data=(node.node_id, row_index, column_index),
                            callback=self._on_lab_cell,
                        )
                        self._lab_cell_ids[(node.node_id, row_index,
                                            column_index)] = cell_id
                        with dpg.item_handler_registry() as registry:
                            dpg.add_item_clicked_handler(
                                button=dpg.mvMouseButton_Left,
                                callback=self._on_lab_cell_click,
                                user_data=(node.node_id, row_index,
                                           column_index),
                            )
                            dpg.add_item_focus_handler(
                                callback=self._on_lab_cell_focus,
                                 user_data=(cell_id, node.node_id,
                                            row_index, column_index),
                                event_type=(dpg.mvEventType_Enter
                                            | dpg.mvEventType_Leave),
                            )
                        dpg.bind_item_handler_registry(cell_id, registry)

    def _ensure_lab_navigation_handlers(self) -> None:
        """Creates one global arrow-key registry for all Lab Data cells."""
        if getattr(self, "_lab_navigation_registry_id", None) is not None:
            return
        with dpg.handler_registry() as registry:
            for key, delta in (
                (dpg.mvKey_Left, (0, -1)),
                (dpg.mvKey_Right, (0, 1)),
                (dpg.mvKey_Up, (-1, 0)),
                (dpg.mvKey_Down, (1, 0)),
            ):
                dpg.add_key_down_handler(
                    key=key, callback=self._on_lab_arrow_key,
                    user_data=delta,
                )
            dpg.add_key_press_handler(
                key=dpg.mvKey_V, callback=self._on_lab_ctrl_v)
        self._lab_navigation_registry_id = registry

    def _on_lab_cell_click(self, sender, app_data, user_data) -> None:
        node_id, row, column = user_data
        self._lab_active_cell = (node_id, row, column)
        if dpg.does_item_exist(sender):
            dpg.bind_item_theme(sender, self._lab_focus_theme())

    def _on_lab_arrow_key(self, sender, app_data, user_data) -> None:
        """Moves focus between measured-data cells with the arrow keys."""
        active = getattr(self, "_lab_active_cell", None)
        if active is None:
            return
        node_id, row, column = active
        current_id = self._lab_cell_ids.get(active)
        if current_id is None or not dpg.does_item_exist(current_id):
            return
        try:
            focused = dpg.get_focused_item()
            lab_ids = set(self._lab_cell_ids.values())
            # DearPyGui can briefly report no focused item while an InputText
            # is handing focus to the next one.  Only reject an arrow from a
            # definitely unrelated widget; this keeps a rapid key sequence
            # moving through the grid instead of stopping after one cell.
            if focused not in (0, None) and focused not in lab_ids:
                return
        except Exception:  # pragma: no cover - defensive for headless DPG
            focused = current_id
        focused_active = next(
            (key for key, value in self._lab_cell_ids.items()
             if value == focused),
            None,
        )
        if focused_active is not None:
            active = focused_active
            node_id, row, column = active
        row_delta, column_delta = user_data
        node = self._state.node_by_id(node_id)
        if node is None:
            return
        rows = self._lab_rows(node, self._lab_columns(node))
        target = (node_id, row + row_delta, column + column_delta)
        if (target[1] < 0 or target[1] >= len(rows)
                or target[2] < 0 or target[2] >= len(self._lab_columns(node))):
            return
        target_id = self._lab_cell_ids.get(target)
        if target_id is None or not dpg.does_item_exist(target_id):
            return
        self._lab_active_cell = target
        dpg.focus_item(target_id)
        dpg.bind_item_theme(target_id, self._lab_focus_theme())

    def _on_lab_ctrl_v(self, sender, app_data, user_data) -> None:
        """Expands a native Ctrl+V directly in the focused Lab Data cell."""
        try:
            if not (dpg.is_key_down(dpg.mvKey_ModCtrl)
                    or dpg.is_key_down(dpg.mvKey_LControl)
                    or dpg.is_key_down(dpg.mvKey_RControl)):
                return
            focused = dpg.get_focused_item()
            target = next(
                (key for key, value in self._lab_cell_ids.items()
                 if value == focused),
                None,
            )
            if target is None:
                return
            node_id, row, column = target
            self._lab_active_cell = target
            # The native InputText paste may already have changed the widget
            # value, but the clipboard remains the canonical range.  If its
            # normal callback has rebuilt the table, ``target`` is no longer
            # present and the method has already returned above.
            self._paste_lab_text(
                node_id, dpg.get_clipboard_text() or "",
                start_row=row, start_column=column)
        except Exception as exc:  # noqa: BLE001 — clipboard/input edge case
            self._set_status(f"Paste failed: {exc}")

    def _rebuild_lab_data_table(self, node_id: str) -> None:
        node = self._state.node_by_id(node_id)
        holder = self._lab_data_holder.get(node_id)
        if node is None or holder is None:
            return
        columns = self._lab_columns(node)
        rows = self._lab_rows(node, columns)
        self._render_lab_data_table(node, holder, columns, rows)
        controls = self._lab_data_controls.get(node_id)
        if controls:
            remove_id, clear_id, count_id = controls
            dpg.configure_item(remove_id, enabled=bool(rows))
            dpg.configure_item(clear_id, enabled=bool(rows))
            dpg.set_value(count_id, f"{len(rows)} point(s)")

    def _on_lab_add_row(self, sender, app_data, user_data) -> None:
        node = self._state.node_by_id(user_data)
        if node is None:
            return
        self._state.add_lab_data_row(user_data, self._lab_columns(node))
        self._rebuild_lab_data_table(user_data)
        self._rebuild_exp_chart_grid(user_data)
        self._schedule_session_autosave()

    def _paste_lab_text(
        self,
        node_id: str,
        text: str,
        *,
        start_row: int | None = None,
        start_column: int | None = None,
    ) -> bool:
        """Applies clipboard text to Lab Data and rebuilds its table."""
        if not str(text).strip():
            self._set_status(
                "Clipboard is empty. Copy an Excel column or table first.")
            return False
        node = self._state.node_by_id(node_id)
        if node is None:
            return False
        columns = self._lab_columns(node)
        active = getattr(self, "_lab_active_cell", None)
        if start_column is None:
            start_column = 0
            if active is not None and active[0] == node_id:
                start_column = max(0, min(active[2], len(columns) - 1))
        if start_row is None:
            start_row = 0
            if active is not None and active[0] == node_id:
                row_count = len(self._lab_rows(node, columns))
                start_row = max(0, min(active[1], row_count))
        rows, had_header = clipboard_service.parse_lab_clipboard(
            str(text), columns, start_column=start_column)
        if not rows:
            self._set_status("Clipboard has no numeric lab data.")
            return False
        self._state.paste_lab_data_rows(node_id, columns, start_row, rows)
        self._lab_active_cell = (node_id, start_row, start_column)
        self._rebuild_lab_data_table(node_id)
        self._rebuild_exp_chart_grid(node_id)
        self._schedule_session_autosave()
        header_note = " with header" if had_header else ""
        self._set_status(
            f"Pasted {len(rows)} lab point(s){header_note} from clipboard.")
        return True

    def _on_lab_paste(self, sender, app_data, user_data) -> None:
        try:
            text = dpg.get_clipboard_text() or ""
        except Exception as exc:  # noqa: BLE001 — clipboard недоступен в headless
            self._set_status(f"Could not read clipboard: {exc}")
            return
        try:
            self._paste_lab_text(user_data, str(text))
        except Exception as exc:  # noqa: BLE001 — callback must not crash the app
            self._set_status(f"Paste failed: {exc}")

    def _on_lab_remove_row(self, sender, app_data, user_data) -> None:
        self._state.remove_lab_data_row(user_data)
        self._rebuild_lab_data_table(user_data)
        self._rebuild_exp_chart_grid(user_data)
        self._schedule_session_autosave()

    def _on_lab_clear(self, sender, app_data, user_data) -> None:
        self._state.clear_lab_data(user_data)
        self._rebuild_lab_data_table(user_data)
        self._rebuild_exp_chart_grid(user_data)
        self._schedule_session_autosave()

    def _on_lab_cell(self, sender, app_data, user_data) -> None:
        node_id, row, column = user_data
        self._lab_active_cell = (node_id, row, column)
        raw = str(app_data).strip()
        multi_value_comma = raw.count(",") >= 2
        if (any(separator in raw for separator in ("\n", "\r", "\t", ";"))
                or ("," in raw and " " in raw) or multi_value_comma):
            try:
                if self._paste_lab_text(node_id, raw, start_row=row,
                                        start_column=column):
                    return
            except Exception as exc:  # noqa: BLE001 — pasted text in one cell
                self._set_status(f"Paste failed: {exc}")
                return
        if not raw:
            value = None
        else:
            try:
                value = float(raw)
                if not math.isfinite(value):
                    raise ValueError("value must be finite")
            except ValueError as exc:
                self._set_status(f"Invalid lab value: {exc}")
                return
        self._state.set_lab_data_value(node_id, row, column, value)
        self._on_lab_cell_focus(sender, None, (sender, node_id, row, column))
        self._rebuild_exp_chart_grid(node_id)
        self._schedule_session_autosave()

    def _lab_focus_theme(self):
        theme_id = getattr(self, "_lab_focus_theme_id", None)
        if theme_id is None:
            with dpg.theme() as theme_id:
                with dpg.theme_component(dpg.mvInputText):
                    dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (255, 244, 180))
                    dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered, (255, 232, 130))
                    dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive, (255, 220, 100))
                    dpg.add_theme_color(dpg.mvThemeCol_InputTextCursor, (30, 30, 30))
            self._lab_focus_theme_id = theme_id
        return theme_id

    def _on_lab_cell_focus(self, sender, app_data, user_data) -> None:
        if isinstance(user_data, tuple) and len(user_data) == 4:
            cell_id, node_id, row, column = user_data
        else:
            cell_id = user_data
            node_id = row = column = None
        if not dpg.does_item_exist(cell_id):
            return
        if dpg.is_item_focused(cell_id):
            if (isinstance(node_id, str) and isinstance(row, int)
                    and isinstance(column, int)):
                self._lab_active_cell = (node_id, row, column)
            dpg.bind_item_theme(cell_id, self._lab_focus_theme())
        else:
            dpg.bind_item_theme(cell_id, 0)

    def _render_exp_table(self, result, parent, cols) -> None:
        allcols = result["columns"]
        idxs = [(c, allcols.index(c)) for c in cols if c in allcols]
        rows = [[self._fmt(row[i]) for _, i in idxs] for row in result["rows"]]
        dpg.add_button(label="Copy table", parent=parent,
                       callback=lambda: self._copy_table(
                           [c for c, _ in idxs], rows, "Experiment results"))
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True, height=-1,
                       freeze_rows=1, freeze_columns=1):
            for c, _ in idxs:
                dpg.add_table_column(label=c)
            for row in rows:
                with dpg.table_row():
                    for value in row:
                        dpg.add_text(value)

    def _render_exp_composition(self, result, parent) -> None:
        stages = result.get("stages", [])
        if not stages:
            dpg.add_text("Per-stage compositions are not available.", parent=parent)
            return
        with dpg.tab_bar(parent=parent):
            for phase in ("liquid", "vapor"):
                with dpg.tab(label=phase.capitalize()) as t:
                    self._render_exp_composition_phase(stages, phase, t)

    def _render_exp_composition_phase(self, stages, phase, parent) -> None:
        names: list = []
        for st in stages:
            if isinstance(st.get(phase), dict):
                names = list(st[phase].keys())
                break
        if not names:
            dpg.add_text(f"No {phase} composition.", parent=parent)
            return
        columns = ["Component"] + [f"{self._g(st.get('pressure'))} bar"
                                   for st in stages]
        rows = []
        for name in names:
            row = [name]
            for st in stages:
                d = st.get(phase)
                row.append(self._fmt(d.get(name) if isinstance(d, dict) else None))
            rows.append(row)
        dpg.add_button(label="Copy table", parent=parent,
                       callback=lambda: self._copy_table(
                           columns, rows, f"{phase.capitalize()} composition"))
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True, height=-1,
                       freeze_rows=1, freeze_columns=1):
            dpg.add_table_column(label="Component")
            for column in columns[1:]:
                dpg.add_table_column(label=column)
            for row in rows:
                with dpg.table_row():
                    for value in row:
                        dpg.add_text(value)

    def _render_exp_chart(self, result, parent, nid) -> None:
        node = self._state.node_by_id(nid)
        with dpg.group(horizontal=True, parent=parent):
            dpg.add_combo(items=result.get("plot_all", []), width=240,
                          label="Add chart for column", user_data=nid,
                          callback=self._on_exp_add_chart)
        holder = dpg.add_group(parent=parent)
        self._exp_chart_holder[nid] = holder
        self._render_exp_chart_grid(result, node, holder)

    def _render_exp_chart_grid(self, result, node, holder) -> None:
        """Рендерит графики карточками по два в строке."""
        if dpg.does_item_exist(holder):
            dpg.delete_item(holder, children_only=True)
        charts = list(result.get("charts", []))
        extra = node.params.get("extra_charts", []) if node else []
        if isinstance(extra, list):
            charts.extend(c for c in extra if c not in charts)
        if not charts:
            dpg.add_text("No plottable columns.", parent=holder)
            return
        lab_data = self._lab_chart_data(node)
        for offset in range(0, len(charts), 2):
            with dpg.group(horizontal=True, parent=holder) as row:
                for col in charts[offset:offset + 2]:
                    with dpg.child_window(width=440, height=270,
                                          border=True, parent=row) as card:
                        self._add_one_chart(result, col, card, lab_data)

    def _rebuild_exp_chart_grid(self, node_id: str) -> None:
        node = self._state.node_by_id(node_id)
        holder = self._exp_chart_holder.get(node_id)
        if (node is None or not isinstance(node.result, dict)
                or holder is None or not dpg.does_item_exist(holder)):
            return
        self._render_exp_chart_grid(node.result, node, holder)

    def _lab_chart_data(self, node):
        if node is None:
            return None
        columns = self._lab_columns(node)
        return {"columns": columns, "rows": self._lab_rows(node, columns),
                "x": "pressure"}

    @staticmethod
    def _exp_x_range(result, lab_result=None) -> tuple[float, float] | None:
        """Диапазон оси X (давления) по данным эксперимента, с малым отступом."""
        xv: list[float] = []
        for source in (result, lab_result):
            if not isinstance(source, dict):
                continue
            cols = source.get("columns", [])
            x = source.get("x")
            if x not in cols:
                continue
            xi = cols.index(x)
            for row in source.get("rows", []):
                if not isinstance(row, list) or xi >= len(row):
                    continue
                value = row[xi]
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    if math.isfinite(float(value)):
                        xv.append(float(value))
        if not xv:
            return None
        lo, hi = min(xv), max(xv)
        margin = (hi - lo) * 0.02 if hi > lo else 1.0
        return lo - margin, hi + margin

    def _add_one_chart(self, result, col, parent, lab_result=None) -> None:
        xs, ys = exp_svc.series_for_plot(result, col)
        lab_x, lab_y = (exp_svc.series_for_plot(lab_result, col)
                        if lab_result else ([], []))
        if not xs and not lab_x:
            return
        xr = self._exp_x_range(result, lab_result)
        with dpg.plot(label=f"{col} vs pressure", height=235, width=-1,
                      parent=parent):
            dpg.add_plot_legend()
            xax = dpg.add_plot_axis(dpg.mvXAxis, label="Pressure, bar")
            yax = dpg.add_plot_axis(dpg.mvYAxis, label=col)
            if xs:
                dpg.add_line_series(xs, ys, label=col, parent=yax)
            if lab_x:
                dpg.add_scatter_series(lab_x, lab_y, label=f"{col} (lab)",
                                       parent=yax)
            if xr is not None:
                dpg.set_axis_limits(xax, xr[0], xr[1])

    def _on_exp_add_chart(self, sender, app_data, user_data) -> None:
        nid = user_data
        col = app_data
        node = self._state.node_by_id(nid)
        if not col or node is None or not isinstance(node.result, dict):
            return
        extra = node.params.get("extra_charts")
        if not isinstance(extra, list):
            extra = []
            node.params["extra_charts"] = extra
        if col not in extra and col not in node.result.get("charts", []):
            extra.append(col)
        self._rebuild_exp_chart_grid(nid)
        self._schedule_session_autosave()
        self._set_status(f"Added chart: {col}")

    def _on_experiment_run(self, sender, app_data, user_data) -> None:
        if self._jobs.busy:
            self._set_status("Another calculation is already running.")
            return
        nid = user_data
        node = self._state.node_by_id(nid)
        composition = self._state.active_composition
        if node is None or composition is None:
            return
        ids = self._exp_input_ids.get(nid, {})
        try:
            pressures = self._parse_floats(dpg.get_value(ids["pressures"]))
            if len(pressures) < 2:
                raise ValueError("need at least 2 pressure stages")
            params = dict(node.params)
            params["pressures"] = pressures
            params["T_c"] = float(dpg.get_value(ids["T_c"]))
            if "P_res" in ids:
                params["P_res"] = float(dpg.get_value(ids["P_res"]))
            if "stage_temps" in ids:
                temps = self._parse_floats(dpg.get_value(ids["stage_temps"]))
                if len(temps) != len(pressures):
                    raise ValueError("stage T count must match pressure count")
                params["stage_temps_c"] = temps
        except (ValueError, KeyError) as exc:
            self._state.set_node_error(nid, f"Invalid input: {exc}")
            return
        self._state.update_node_params(nid, params, notify=False)
        kind = node.params["kind"]

        ref = self._state.node_ref(nid)
        if ref is None:
            return
        self._jobs.start(ref, kind,
                         lambda: exp_svc.run_experiment(composition, kind, params))
        self._state.set_node_running(ref)
        self._set_status(f"{kind.upper()} running...")
        self._arm_flash_poll()

    @staticmethod
    def _parse_floats(text: str) -> list[float]:
        """Парсит числа из строки (разделители — запятая/пробел/перенос строки)."""
        out: list[float] = []
        for chunk in text.replace("\n", ",").replace(" ", ",").split(","):
            chunk = chunk.strip()
            if chunk:
                out.append(float(chunk))
        return out
