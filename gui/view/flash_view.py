"""Вкладки Flash/Compare и callbacks фонового flash-расчёта."""

import dearpygui.dearpygui as dpg

from gui.app_state import NodeStatus
from gui.services import flash_service
from gui.view.contracts import ContextBoundView
from gui.view.read_only_table import render_readonly_table


class FlashViewMixin(ContextBoundView):
    """Рендер Flash/Compare и приём завершённых GUI-задач."""

    _flash_input_ids: dict[str, tuple[int, int]]

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
        dpg.add_button(label="Copy table", parent=parent,
                       callback=lambda: self._copy_table(
                           ["Property", "Vapor", "Liquid"], rows,
                           "Flash results"))
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
        dpg.add_button(label="Copy table", parent=parent,
                       callback=lambda: self._copy_table(
                           ["Component", "Vapor yi", "Liquid xi", "K = yi/xi"],
                           rows, "Flash composition"))
        render_readonly_table(
            parent, ["Component", "Vapor yi", "Liquid xi", "K = yi/xi"], rows)

    # ==================================================================
    #  Вкладка Compare (таблица сравнения + плитка панелей)
    # ==================================================================

    def _render_compare_tab(self, parent, node) -> None:
        members = [self._state.node_by_id(m) for m in node.params.get("members", [])]
        members = [m for m in members if m is not None and m.result is not None]
        if len(members) < 2:
            dpg.add_text("Select at least 2 computed flash runs to compare "
                         "(Ctrl+click runs in the tree, or use the right-click menu).",
                         parent=parent)
            return
        headers = [self._compare_col_label(m) for m in members]
        with dpg.tab_bar(parent=parent):
            with dpg.tab(label="Vapor") as t:
                self._compare_phase_table(members, headers, "vapor", t)
            with dpg.tab(label="Liquid") as t:
                self._compare_phase_table(members, headers, "liquid", t)
            with dpg.tab(label="K-values") as t:
                self._compare_k_table(members, headers, t)
            with dpg.tab(label="Panels") as t:
                self._compare_panels(members, headers, t)

    def _compare_col_label(self, node) -> str:
        label = node.params.get("label")
        if label:
            return label
        return f"{self._g(node.params.get('P'))}/{self._g(node.params.get('T'))}"

    def _compare_phase_table(self, members, headers, phase, parent) -> None:
        rows = [["Phase mole fraction"] + [
            self._fmt(getattr(member.result, phase).mole_fraction)
            for member in members]]
        rows.extend([
            [label] + [self._fmt(getattr(member.result, phase).properties.get(key))
                       for member in members]
            for key, label in flash_service.PHASE_PROPERTY_ROWS
        ])
        dpg.add_button(label="Copy table", parent=parent,
                       callback=lambda: self._copy_table(
                           ["Property"] + headers, rows, "Comparison"))
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
        dpg.add_button(label="Copy table", parent=parent,
                       callback=lambda: self._copy_table(
                           ["Component (K=yi/xi)"] + headers, rows,
                           "K-value comparison"))
        render_readonly_table(parent, ["Component (K=yi/xi)"] + headers, rows)

    def _compare_panels(self, members, headers, parent) -> None:
        with dpg.group(horizontal=True, parent=parent):
            for m, h in zip(members, headers):
                with dpg.child_window(width=380, height=-1, border=True) as cw:
                    dpg.add_text(h)
                    self._render_flash_result(m.result, cw)

    def _flash_pt(self, nid: str) -> tuple[float, float]:
        """Текущие значения полей P/T вкладки флэша по её id узла."""
        pid, tid = self._flash_input_ids[nid]
        return dpg.get_value(pid), dpg.get_value(tid)

    # --- фоновая задача флэша (поток + опрос по кадрам) ------------------

    def _on_flash_param(self, sender, app_data, user_data) -> None:
        nid = user_data
        p, t = self._flash_pt(nid)
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
        p, t = self._flash_pt(nid)
        self._state.set_flash_params(nid, p, t)

        ref = self._state.node_ref(nid)
        if ref is None:
            return
        self._jobs.start(ref, "flash",
                         lambda: flash_service.run_flash(composition, p, t))
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
