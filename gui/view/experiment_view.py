"""Вкладки PVT-экспериментов, таблицы и графики результатов."""

import math

import dearpygui.dearpygui as dpg

from gui.app_state import NodeStatus
from gui.services import clipboard_service, input_validation_service
from gui.services import experiment_service as exp_svc
from gui.services import lab_data_service as lab_svc
from gui.view.contracts import ContextBoundView
from gui.view.read_only_table import render_readonly_table


class ExperimentViewMixin(ContextBoundView):
    """Рендер и callbacks CCE/DLE/Separator."""

    _exp_input_ids: dict[str, dict[str, int]]
    _exp_validation_message_ids: dict[str, int]
    _exp_chart_holder: dict[str, int]
    _lab_data_holder: dict[str, int]
    _lab_data_controls: dict[str, tuple[int, int, int]]
    _lab_source_choices: dict[str, dict[str, str | None]]
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
            default_value=", ".join(self._g(x) for x in p.get("pressures", [])),
            user_data=nid, callback=self._on_experiment_input_changed)
        with dpg.group(horizontal=True, parent=parent):
            ids["T_c"] = dpg.add_input_float(label="T, C", width=140, step=0,
                                             default_value=float(p.get("T_c", 100.0)),
                                             user_data=nid,
                                             callback=self._on_experiment_input_changed)
            if meta.get("needs_p_res"):
                ids["P_res"] = dpg.add_input_float(
                    label="P res, bar", width=140, step=0,
                    default_value=float(p.get("P_res", 400.0)), user_data=nid,
                    callback=self._on_experiment_input_changed)
        if meta.get("needs_stage_temps"):
            ids["stage_temps"] = dpg.add_input_text(
                label="Stage T, C (comma-separated, same count)", width=460,
                parent=parent,
                default_value=", ".join(self._g(x) for x in p.get("stage_temps_c", [])),
                user_data=nid, callback=self._on_experiment_input_changed)
        self._exp_input_ids[nid] = ids
        self._exp_validation_message_ids[nid] = dpg.add_text(
            "", parent=parent, show=False, wrap=620)

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
        data = self._effective_lab_data(node)
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

    def _effective_lab_data(self, node):
        model = self._state.active_model
        if model is None:
            return node.params.get("lab_data")
        return lab_svc.effective_lab_data(
            self._state.db_path, model.project_id, model.model_id, node.params,
        )

    def _linked_lab_dataset(self, node):
        ref = node.params.get("lab_data_ref")
        model = self._state.active_model
        if not isinstance(ref, str) or not ref or model is None:
            return None
        try:
            return lab_svc.get_dataset(self._state.db_path, model.project_id, ref,
                                       model_id=model.model_id)
        except lab_svc.LabDataStoreError:
            return None

    @staticmethod
    def _lab_condition_value(raw) -> float | None:
        if isinstance(raw, bool) or not isinstance(raw, (int, float)):
            return None
        value = float(raw)
        return value if math.isfinite(value) else None

    def _render_linked_lab_context(self, node, dataset, parent) -> None:
        """Shows non-blocking context for a selected project Lab Data set."""
        conditions = dataset.get("conditions")
        conditions = conditions if isinstance(conditions, dict) else {}
        displayed: list[str] = []
        mismatches: list[str] = []
        for param, condition_key, label, unit in (
            ("T_c", "T_c", "T res", "C"),
            ("P_res", "P_res", "P res", "bar"),
        ):
            lab_value = self._lab_condition_value(conditions.get(condition_key))
            if lab_value is None:
                continue
            displayed.append(f"{label} = {self._g(lab_value)} {unit}")
            calculation_value = self._lab_condition_value(node.params.get(param))
            if (calculation_value is not None
                    and not math.isclose(lab_value, calculation_value,
                                         rel_tol=0.0, abs_tol=1e-6)):
                mismatches.append(
                    f"{label}: Lab Data {self._g(lab_value)} {unit}; "
                    f"calculation {self._g(calculation_value)} {unit}."
                )
        if displayed:
            dpg.add_text("Dataset conditions: " + "; ".join(displayed) + ".",
                         parent=parent)
        else:
            note = dpg.add_text("Dataset conditions are not specified.", parent=parent)
            dpg.bind_item_theme(note, self._theme_stale())
        for message in mismatches:
            warning = dpg.add_text("Warning — " + message, parent=parent, wrap=620)
            dpg.bind_item_theme(warning, self._theme_stale())

        model = self._state.active_model
        if model is None:
            return
        references = lab_svc.dataset_references(
            self._state.db_path, model.project_id, dataset["dataset_id"])
        references.discard((model.model_id, node.node_id))
        dpg.add_text(
            f"Used by {len(references)} other saved calculation(s) in this project.",
            parent=parent,
        )

    def _render_lab_data(self, node, parent) -> None:
        """Сворачиваемый ввод измерений, привязанный к одному Experiment-узлу."""
        columns = self._lab_columns(node)
        rows = self._lab_rows(node, columns)
        linked = self._linked_lab_dataset(node)
        with dpg.collapsing_header(label="Lab Data (measured)",
                                  default_open=True, parent=parent) as header:
            model = self._state.active_model
            choices: dict[str, str | None] = {"No Lab Data": None}
            if model is not None:
                try:
                    datasets = lab_svc.list_datasets(
                        self._state.db_path, model.project_id,
                        experiment_kind=node.params.get("kind"))
                except lab_svc.LabDataStoreError as exc:
                    datasets = []
                    dpg.add_text(str(exc), parent=header)
                for dataset in datasets:
                    choices[dataset["title"]] = dataset["dataset_id"]
            self._lab_source_choices[node.node_id] = choices
            selected = next((label for label, value in choices.items()
                             if linked and value == linked["dataset_id"]),
                            "No Lab Data")
            dpg.add_combo(label="Data source", items=list(choices), width=360,
                          default_value=selected, user_data=node.node_id,
                          callback=self._on_lab_source_changed, parent=header)
            if linked is not None:
                dpg.add_text(
                    f"{linked['scope'].title()} dataset: {linked['title']} "
                    "(read only).",
                    parent=header,
                )
                self._render_linked_lab_context(node, linked, header)
                dpg.add_button(label="Edit selected Lab Data",
                               user_data=node.node_id,
                               callback=self._on_lab_open_selected,
                               parent=header)
            else:
                message = ("Choose a Lab Data source in the list above."
                           if not rows else
                           "Legacy local Lab Data is read only; choose a catalog source.")
                dpg.add_text(message, parent=header, wrap=620)
            dpg.add_text(f"{len(rows)} point(s)", parent=header)
            holder = dpg.add_group(parent=header)
            self._lab_data_holder[node.node_id] = holder
            self._render_lab_data_table(node, holder, columns, rows,
                                        readonly=True)

    def _render_lab_data_table(self, node, holder, columns, rows, *, readonly=False) -> None:
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
                        if readonly:
                            dpg.add_text("" if value is None else self._g(value))
                            continue
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
        self._render_lab_data_table(node, holder, columns, rows,
                                    readonly=self._linked_lab_dataset(node) is not None)
        controls = self._lab_data_controls.get(node_id)
        if controls:
            remove_id, clear_id, count_id = controls
            dpg.configure_item(remove_id, enabled=bool(rows))
            dpg.configure_item(clear_id, enabled=bool(rows))
            dpg.set_value(count_id, f"{len(rows)} point(s)")

    def _on_lab_source_changed(self, sender, app_data, user_data) -> None:
        node_id = str(user_data)
        choices = self._lab_source_choices.get(node_id, {})
        self._state.set_lab_data_ref(node_id, choices.get(str(app_data)))
        self._render_workspace()
        self._schedule_session_autosave()

    def _on_lab_open_selected(self, sender, app_data, user_data) -> None:
        node = self._state.node_by_id(str(user_data))
        if node is None:
            return
        dataset = self._linked_lab_dataset(node)
        if dataset is not None:
            self._open_lab_dataset_editor(dataset)

    def _on_lab_make_local_copy(self, sender, app_data, user_data) -> None:
        node_id = str(user_data)
        node = self._state.node_by_id(node_id)
        if node is None:
            return
        data = self._effective_lab_data(node)
        if not isinstance(data, dict):
            return
        columns = data.get("columns")
        rows = data.get("rows")
        if not isinstance(columns, list) or not isinstance(rows, list):
            return
        self._state.copy_lab_data_to_local(node_id, list(columns), list(rows))
        self._render_workspace()
        self._schedule_session_autosave()
        self._set_status("Lab Data copied to this run for editing.")

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
        if not dpg.does_item_exist(sender):
            return
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
                # Respect the decimal-comma locale used by Excel/Russian
                # spreadsheets for ordinary one-cell manual edits.
                scalar = raw.replace(",", ".") if raw.count(",") == 1 else raw
                value = float(scalar)
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
        self._add_table_copy_controls(
            parent, [c for c, _ in idxs], rows, "Experiment results")
        render_readonly_table(parent, [c for c, _ in idxs], rows)

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
        self._add_table_copy_controls(
            parent, columns, rows, f"{phase.capitalize()} composition")
        render_readonly_table(parent, columns, rows)

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
        cards_per_row = self._chart_grid_columns()
        for offset in range(0, len(charts), cards_per_row):
            with dpg.group(horizontal=(cards_per_row > 1), parent=holder) as row:
                for col in charts[offset:offset + cards_per_row]:
                    with dpg.child_window(
                            width=self._chart_card_width(cards_per_row),
                            height=self._chart_card_height(),
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
        with dpg.plot(label=f"{col} vs pressure", height=self._chart_plot_height(),
                      width=-1,
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
        if not col:
            return
        if not self._state.add_experiment_chart(nid, col):
            return
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
        params = self._validated_experiment_params(nid, node.params.get("kind", "cce"))
        if params is None:
            self._set_status("Correct the highlighted experiment inputs before running.")
            return
        self._state.update_node_params(nid, params, notify=False)
        kind = node.params["kind"]

        ref = self._state.node_ref(nid)
        if ref is None:
            return
        self._jobs.start(
            ref, kind,
            lambda token, progress: exp_svc.run_experiment(
                composition, kind, params,
                cancellation_token=token,
                progress_callback=progress,
            ),
        )
        self._state.set_node_running(ref)
        self._set_status(f"{kind.upper()} running...")
        self._arm_flash_poll()

    def _validated_experiment_params(self, nid: str, kind: str) -> dict | None:
        ids = self._exp_input_ids.get(nid, {})
        if not ids:
            return None
        parsed, errors = input_validation_service.validate_experiment_inputs(
            kind,
            dpg.get_value(ids["pressures"]),
            dpg.get_value(ids["T_c"]),
            dpg.get_value(ids["P_res"]) if "P_res" in ids else None,
            dpg.get_value(ids["stage_temps"]) if "stage_temps" in ids else None,
        )
        if not self._show_input_validation(
                ids, errors, self._exp_validation_message_ids.get(nid)):
            return None
        node = self._state.node_by_id(nid)
        return {**(node.params if node is not None else {}), **parsed}

    def _on_experiment_input_changed(self, sender, app_data, user_data) -> None:
        node = self._state.node_by_id(str(user_data))
        if node is not None:
            self._validated_experiment_params(str(user_data), node.params.get("kind", "cce"))

    @staticmethod
    def _parse_floats(text: str) -> list[float]:
        """Парсит числа из строки (разделители — запятая/пробел/перенос строки)."""
        out: list[float] = []
        for chunk in text.replace("\n", ",").replace(" ", ",").split(","):
            chunk = chunk.strip()
            if chunk:
                out.append(float(chunk))
        return out
