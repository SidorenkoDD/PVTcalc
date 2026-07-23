"""Вкладки Flash/Compare и callbacks фонового flash-расчёта."""

import math

import dearpygui.dearpygui as dpg

from gui.app_state import NodeKind, NodeRef, NodeStatus
from gui.services import comparison_service, flash_service, input_validation_service
from gui.services import lab_data_service as lab_svc
from gui.view.contracts import ContextBoundView
from gui.view.read_only_table import render_readonly_table


class FlashViewMixin(ContextBoundView):
    """Рендер Flash/Compare и приём завершённых GUI-задач."""

    _flash_input_ids: dict[str, tuple[int, int]]
    _flash_validation_message_ids: dict[str, int]

    def _render_flash_tab(self, parent, node) -> None:
        nid = node.node_id
        params = node.params
        running = node.status is NodeStatus.RUNNING

        with dpg.group(horizontal=True, parent=parent):
            pid = dpg.add_input_float(label="P, bar", width=150, step=0,
                                      default_value=float(params.get("P", 100.0)),
                                      user_data=nid, callback=self._on_flash_param)
            tid = dpg.add_input_float(label="T, C", width=150, step=0,
                                      default_value=float(params.get("T", 20.0)),
                                      user_data=nid, callback=self._on_flash_param)
            self._flash_input_ids[nid] = (pid, tid)
            if running:
                dpg.add_loading_indicator(style=1, radius=2.0)
                dpg.add_button(label="Cancel", user_data=nid,
                               callback=self._on_flash_cancel)
            else:
                dpg.add_button(label="Run flash", user_data=nid,
                               callback=self._on_flash_run)
        self._flash_validation_message_ids[nid] = dpg.add_text(
            "", parent=parent, show=False, wrap=620)

        if running:
            dpg.add_text("Running flash...", parent=parent)
        if node.status is NodeStatus.STALE and node.error:
            err = dpg.add_text(f"Flash error: {node.error}", parent=parent)
            dpg.bind_item_theme(err, self._theme_stale())
        elif node.status is NodeStatus.STALE and node.result is not None:
            st = dpg.add_text("Result is stale (composition changed) - re-run.",
                              parent=parent)
            dpg.bind_item_theme(st, self._theme_stale())

        if node.result is not None:
            self._render_flash_result(node.result, parent)

    def _render_flash_result(self, result, parent) -> None:
        dpg.add_separator(parent=parent)
        kind = "Two-phase" if result.is_two_phase else "Single-phase"
        dpg.add_text(f"{kind}    P = {result.pressure:.3f} bar    "
                     f"T = {result.temperature:.2f} K", parent=parent)
        if result.diagnostics.warnings:
            warning = dpg.add_text(
                "Result obtained, but requires verification:\n" + "\n".join(
                    f"- {item.message}" for item in result.diagnostics.warnings
                ),
                wrap=760,
                parent=parent,
            )
            dpg.bind_item_theme(warning, self._theme_stale())
        # незакрываемые под-вкладки результата: Main + Composition
        with dpg.tab_bar(parent=parent):
            with dpg.tab(label="Main") as t_main:
                self._render_flash_main(result, t_main)
            with dpg.tab(label="Composition & K") as t_comp:
                self._render_flash_composition(result, t_comp)

    def _render_flash_main(self, result, parent) -> None:
        """Основные свойства фаз (доли, M, объём, плотность, Z, вязкость)."""
        rows = [["Phase mole fraction", self._fmt(result.vapor.mole_fraction),
                 self._fmt(result.liquid.mole_fraction)]]
        rows.extend([
            [label, self._fmt(result.vapor.properties.get(key)),
             self._fmt(result.liquid.properties.get(key))]
            for key, label in flash_service.PHASE_PROPERTY_ROWS
        ])
        self._add_table_copy_controls(
            parent, ["Property", "Vapor", "Liquid"], rows, "Flash results")
        render_readonly_table(parent, ["Property", "Vapor", "Liquid"], rows,
                              scroll_x=False, freeze_columns=0)

    def _render_flash_composition(self, result, parent) -> None:
        """Состав каждой фазы (yi/xi) и константы равновесия K = yi/xi."""
        yi = result.vapor.composition if isinstance(result.vapor.composition, dict) else None
        xi = result.liquid.composition if isinstance(result.liquid.composition, dict) else None
        if not yi and not xi:
            dpg.add_text("Phase compositions are not available for this result.",
                         parent=parent)
            return

        phase_composition = xi if xi else yi
        if phase_composition is None:
            return
        names = list(phase_composition.keys())
        rows = []
        for name in names:
            y = yi.get(name) if yi else None
            x = xi.get(name) if xi else None
            k = (y / x) if (y is not None and x not in (None, 0.0)) else None
            rows.append([name, self._fmt(y), self._fmt(x), self._fmt(k)])
        self._add_table_copy_controls(
            parent, ["Component", "Vapor yi", "Liquid xi", "K = yi/xi"], rows,
            "Flash composition")
        render_readonly_table(
            parent, ["Component", "Vapor yi", "Liquid xi", "K = yi/xi"], rows)

    # ==================================================================
    #  Вкладка Compare (таблица сравнения + плитка панелей)
    # ==================================================================

    def _resolve_compare_members(self, compare_node):
        """Разрешает новые межмодельные refs и старые локальные members."""
        refs = []
        raw_refs = compare_node.params.get("member_refs")
        if isinstance(raw_refs, list):
            refs = [ref for value in raw_refs
                    if (ref := NodeRef.from_dict(value)) is not None]
        if not refs:
            refs = [ref for node_id in compare_node.params.get("members", [])
                    if (ref := self._state.node_ref(str(node_id))) is not None]
        return [(ref, member) for ref in refs
                if (member := self._state.node_by_ref(ref)) is not None]

    def _render_compare_tab(self, parent, node) -> None:
        resolved = self._resolve_compare_members(node)
        computed = [(ref, member) for ref, member in resolved
                    if member.result is not None]
        if len(computed) < 2:
            dpg.add_text("Select at least 2 computed runs to compare "
                         "(Ctrl+click runs in the tree, or use the right-click menu).",
                         parent=parent)
            return
        if all(member.kind is NodeKind.EXPERIMENT for _ref, member in computed):
            self._render_experiment_compare(computed, parent)
            return
        if not all(member.kind is NodeKind.FLASH for _ref, member in computed):
            dpg.add_text("Selected results have incompatible types.", parent=parent)
            return
        members = [member for _ref, member in computed]
        headers = [self._compare_col_label(ref, member) for ref, member in computed]
        with dpg.tab_bar(parent=parent):
            with dpg.tab(label="Vapor") as t:
                self._compare_phase_table(members, headers, "vapor", t)
            with dpg.tab(label="Liquid") as t:
                self._compare_phase_table(members, headers, "liquid", t)
            with dpg.tab(label="K-values") as t:
                self._compare_k_table(members, headers, t)
            with dpg.tab(label="Panels") as t:
                self._compare_panels(members, headers, t)

    def _render_experiment_compare(self, members, parent) -> None:
        descriptors = []
        for ref, member in members:
            label = self._exp_compare_label(ref, member)
            model = self._state.models.get(ref.model_id)
            raw_lab = lab_svc.effective_lab_data(
                self._state.db_path,
                model.project_id if model is not None else "",
                ref.model_id,
                member.params,
            )
            lab_data = None
            if isinstance(raw_lab, dict):
                lab_data = {
                    "columns": raw_lab.get("columns", []),
                    "rows": raw_lab.get("rows", []),
                    "x": "pressure",
                }
            descriptors.append({
                "id": f"{ref.model_id}/{ref.variant_id}/{ref.node_id}",
                "label": label,
                "kind": member.params.get("kind"),
                "params": dict(member.params),
                "result": member.result,
                "lab_data": lab_data,
                "stale": member.status is NodeStatus.STALE,
            })
        try:
            comparison = comparison_service.build_experiment_comparison(descriptors)
        except comparison_service.IncompatibleComparisonError as exc:
            dpg.add_text(str(exc), parent=parent)
            return

        dpg.add_text(
            f"{str(comparison['kind']).upper()} comparison; "
            f"reference: {comparison['reference']}",
            parent=parent,
        )
        for message in comparison["warnings"]:
            warning = dpg.add_text(f"Warning: {message}", wrap=760, parent=parent)
            dpg.bind_item_theme(warning, self._theme_stale())
        with dpg.tab_bar(parent=parent):
            with dpg.tab(label="Curves") as curves:
                self._render_experiment_compare_charts(comparison, curves)
            with dpg.tab(label="Deviations") as deviations:
                self._render_experiment_deviations(comparison, deviations)

    def _exp_compare_label(self, ref: NodeRef, node) -> str:
        model = self._state.models.get(ref.model_id)
        model_label = model.title if model is not None else ref.model_id
        seq = node.node_id.split("_")[-1]
        run_label = (node.params.get("label")
                     or f"{str(node.params.get('kind', 'experiment')).upper()} #{seq}")
        return f"{model_label} · {run_label}"

    def _render_experiment_compare_charts(self, comparison, parent) -> None:
        columns = comparison["columns"]
        cards_per_row = self._chart_grid_columns()
        for offset in range(0, len(columns), cards_per_row):
            with dpg.group(horizontal=(cards_per_row > 1), parent=parent) as row:
                for column in columns[offset:offset + cards_per_row]:
                    with dpg.child_window(
                            width=self._chart_card_width(cards_per_row),
                            height=self._chart_card_height(), border=True,
                                          parent=row) as card:
                        self._add_experiment_compare_chart(
                            column, comparison["series"][column], card,
                        )

    def _add_experiment_compare_chart(self, column, series, parent) -> None:
        all_x = [float(value) for item in series
                 for values in (item["x"], item["lab_x"])
                 for value in values if isinstance(value, (int, float))
                 and not isinstance(value, bool) and math.isfinite(float(value))]
        with dpg.plot(label=f"{column} vs pressure",
                      height=self._chart_plot_height(), width=-1,
                      parent=parent):
            dpg.add_plot_legend()
            x_axis = dpg.add_plot_axis(dpg.mvXAxis, label="Pressure, bar")
            y_axis = dpg.add_plot_axis(dpg.mvYAxis, label=column)
            for item in series:
                if item["x"]:
                    dpg.add_line_series(item["x"], item["y"],
                                        label=item["label"], parent=y_axis)
                if item["lab_x"]:
                    dpg.add_scatter_series(
                        item["lab_x"], item["lab_y"],
                        label=f"{item['label']} (lab)", parent=y_axis,
                    )
            if all_x:
                lo, hi = min(all_x), max(all_x)
                margin = (hi - lo) * 0.02 if hi > lo else 1.0
                dpg.set_axis_limits(x_axis, lo - margin, hi + margin)

    def _render_experiment_deviations(self, comparison, parent) -> None:
        with dpg.tab_bar(parent=parent):
            for column in comparison["columns"]:
                with dpg.tab(label=column) as tab:
                    rows = [
                        [
                            self._fmt(row["pressure"]), row["member"],
                            self._fmt(row["reference"]), self._fmt(row["candidate"]),
                            self._fmt(row["absolute"]),
                            self._fmt(row["relative_percent"]),
                        ]
                        for row in comparison["deviations"][column]
                    ]
                    headers = [
                        "Pressure, bar", "Compared run", "Reference",
                        "Compared", "|difference|", "Difference, %",
                    ]
                    if not rows:
                        dpg.add_text(
                            "No matching pressure points for this property.",
                            parent=tab,
                        )
                        continue
                    self._add_table_copy_controls(
                        tab, headers, rows, f"{column} deviations")
                    render_readonly_table(tab, headers, rows)

    def _compare_col_label(self, ref: NodeRef, node) -> str:
        model = self._state.models.get(ref.model_id)
        model_label = model.title if model is not None else ref.model_id
        run_label = (node.params.get("label")
                     or f"{self._g(node.params.get('P'))}/{self._g(node.params.get('T'))}")
        return f"{model_label} · {run_label}"

    def _compare_phase_table(self, members, headers, phase, parent) -> None:
        rows = [["Phase mole fraction"] + [
            self._fmt(getattr(member.result, phase).mole_fraction)
            for member in members]]
        rows.extend([
            [label] + [self._fmt(getattr(member.result, phase).properties.get(key))
                       for member in members]
            for key, label in flash_service.PHASE_PROPERTY_ROWS
        ])
        self._add_table_copy_controls(
            parent, ["Property"] + headers, rows, "Comparison")
        render_readonly_table(parent, ["Property"] + headers, rows)

    def _compare_k_table(self, members, headers, parent) -> None:
        # набор компонентов — из первого участника с известным составом жидкости
        names: list = []
        for m in members:
            xi = m.result.liquid.composition
            if isinstance(xi, dict):
                names = list(xi.keys())
                break
        if not names:
            dpg.add_text("Phase compositions are not available for these results.",
                         parent=parent)
            return
        rows = []
        for name in names:
            row = [name]
            for member in members:
                yi = member.result.vapor.composition
                xi = member.result.liquid.composition
                y = yi.get(name) if isinstance(yi, dict) else None
                x = xi.get(name) if isinstance(xi, dict) else None
                row.append(self._fmt((y / x) if (y is not None
                                  and x not in (None, 0.0)) else None))
            rows.append(row)
        self._add_table_copy_controls(
            parent, ["Component (K=yi/xi)"] + headers, rows,
            "K-value comparison")
        render_readonly_table(parent, ["Component (K=yi/xi)"] + headers, rows)

    def _compare_panels(self, members, headers, parent) -> None:
        cards_per_row = self._chart_grid_columns()
        for offset in range(0, len(members), cards_per_row):
            with dpg.group(horizontal=(cards_per_row > 1), parent=parent):
                for m, h in zip(members[offset:offset + cards_per_row],
                                headers[offset:offset + cards_per_row]):
                    with dpg.child_window(width=self._chart_card_width(cards_per_row),
                                          height=-1, border=True) as cw:
                        dpg.add_text(h)
                        self._render_flash_result(m.result, cw)

    def _flash_pt(self, nid: str) -> tuple[float, float]:
        """Текущие значения полей P/T вкладки флэша по её id узла."""
        pid, tid = self._flash_input_ids[nid]
        return dpg.get_value(pid), dpg.get_value(tid)

    def _validate_flash_form(self, nid: str) -> tuple[float, float] | None:
        p, t = self._flash_pt(nid)
        pid, tid = self._flash_input_ids[nid]
        errors = input_validation_service.validate_flash_inputs(p, t)
        if not self._show_input_validation(
                {"P": pid, "T": tid}, errors,
                self._flash_validation_message_ids.get(nid)):
            return None
        return p, t

    # --- фоновая задача флэша (поток + опрос по кадрам) ------------------

    def _on_flash_param(self, sender, app_data, user_data) -> None:
        nid = user_data
        values = self._validate_flash_form(nid)
        if values is None:
            return
        p, t = values
        self._state.set_flash_params(nid, p, t)

    def _on_flash_run(self, sender, app_data, user_data) -> None:
        if self._jobs.busy:
            self._set_status("Another flash is already running.")
            return
        nid = user_data
        composition = self._state.active_composition
        node = self._state.node_by_id(nid)
        if composition is None or node is None:
            return
        values = self._validate_flash_form(nid)
        if values is None:
            self._set_status("Correct the highlighted flash inputs before running.")
            return
        p, t = values
        self._state.set_flash_params(nid, p, t)

        ref = self._state.node_ref(nid)
        if ref is None:
            return
        self._jobs.start(
            ref, "flash",
            lambda token, progress: flash_service.run_flash(
                composition, p, t,
                cancellation_token=token,
                progress_callback=progress,
            ),
        )
        self._state.set_flash_running(ref)
        self._set_status(f"Flash running at P={p:.3f} bar, T={t:.2f} C...")
        self._arm_flash_poll()

    def _arm_flash_poll(self) -> None:
        dpg.set_frame_callback(dpg.get_frame_count() + 1, self._poll_flash_job)

    def _poll_flash_job(self, sender=None, app_data=None) -> None:
        job = self._jobs.active
        if job is None:
            return
        if not job.done.is_set():
            if job.progress_message:
                self._set_status(
                    f"{job.label.capitalize()}: {job.progress_message} "
                    f"({job.progress:.0%})"
                )
            self._arm_flash_poll()
            return
        job = self._jobs.take_finished()
        if job is None:
            return
        if job.cancelled:
            self._state.reset_node(job.node_ref)
            self._set_status(f"{job.label.capitalize()} cancelled.")
        elif job.error:
            self._state.set_node_error(job.node_ref, job.error)
            self._set_status(f"{job.label.capitalize()} failed.")
        else:
            self._state.set_node_result(job.node_ref, job.result)
            self._set_status(f"{job.label.capitalize()} done.")

    def _on_flash_cancel(self, sender, app_data, user_data) -> None:
        if self._jobs.cancel():
            self._set_status("Cancelling flash (result will be discarded)...")
