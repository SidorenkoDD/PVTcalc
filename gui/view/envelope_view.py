"""Вкладка фазовой огибающей, её график и диалог параметров."""

import dearpygui.dearpygui as dpg

from gui.app_state import NodeStatus
from gui.services import phase_envelope_service as pe_svc
from gui.view.contracts import ContextBoundView
from gui.view.read_only_table import render_readonly_table

_ENV_METHOD_SSM = "SSM (curve)"
_ENV_METHOD_GRID = "Grid (stability scan)"


def _env_method_from_label(label: str) -> str:
    """Подпись метода из combo -> ключ метода."""
    return "grid" if str(label).startswith("Grid") else "ssm"


class EnvelopeViewMixin(ContextBoundView):
    """Рендер и callbacks Phase Envelope."""

    _env_dialog_win: int | None
    _env_dialog_ids: dict[str, int]
    _env_dialog_node: str | None

    def _render_envelope_tab(self, parent, node) -> None:
        nid = node.node_id
        p = node.params
        running = node.status is NodeStatus.RUNNING

        # параметры — read-only сводка (правятся во всплывающем окне)
        if p.get("method") == "grid":
            dpg.add_text("Phase envelope (P-T): thermodynamic stability scan on "
                         "a P-T grid.", parent=parent)
            summary = (f"Grid: T = [0 .. {self._g(p.get('grid_t_max_c'))}] C, "
                       f"P = [1 .. {self._g(p.get('grid_p_max_bar'))}] bar, "
                       f"{self._g(p.get('grid_t_points'))}x{self._g(p.get('grid_p_points'))} points")
        else:
            dpg.add_text("Phase envelope (P-T): saturation curve (SSM) + reservoir "
                         "saturation pressure.", parent=parent)
            summary = (f"SSM: T = [{self._g(p.get('t_min_c'))} .. {self._g(p.get('t_max_c'))}] C, "
                       f"step {self._g(p.get('t_step_c'))} C    |    "
                       f"P max = {self._g(p.get('p_max_bar'))} bar")
        dpg.add_text(summary, parent=parent)
        tres = p.get("T_res_c")
        if tres is not None:
            dpg.add_text(f"Reservoir T = {self._g(tres)} C (Psat marker)", parent=parent)

        with dpg.group(horizontal=True, parent=parent):
            if running:
                dpg.add_loading_indicator(style=1, radius=2.0)
                dpg.add_text("Parallel calculation cannot be cancelled.")
            else:
                dpg.add_button(label="Edit parameters & run...", user_data=nid,
                               callback=self._on_envelope_edit)
        if running:
            dpg.add_text("Computing envelope in parallel (this may take up to a minute)...",
                         parent=parent)
        if node.status is NodeStatus.STALE and node.error:
            err = dpg.add_text(f"Error: {node.error}", parent=parent)
            dpg.bind_item_theme(err, self._theme_stale())
        elif node.status is NodeStatus.STALE and node.result is not None:
            st = dpg.add_text("Result is stale (composition changed) - re-run.",
                              parent=parent)
            dpg.bind_item_theme(st, self._theme_stale())

        if node.result is not None:
            dpg.add_separator(parent=parent)
            with dpg.tab_bar(parent=parent):
                with dpg.tab(label="Chart") as t:
                    self._render_envelope_chart(node.result, t)
                with dpg.tab(label="Data") as t:
                    self._render_envelope_table(node.result, t)

    def _render_envelope_chart(self, result, parent) -> None:
        res = result.get("reservoir")
        if result.get("method") == "grid":
            grid = result.get("grid", {})
            uns, sta = grid.get("unstable", {}), grid.get("stable", {})
            all_p = [v for v in (uns.get("P", []) + sta.get("P", [])) if v is not None]
        else:
            env = result.get("envelope", {})
            all_p = [v for v in env.get("P", []) if v is not None]
        if res and res.get("P_sat") is not None:
            all_p.append(res["P_sat"])
        if not all_p:
            dpg.add_text("No envelope data - try a wider range or higher P max.",
                         parent=parent)
            return
        p_hi = max(all_p) * 1.05

        # фиксированная высота — надёжнее, чем height=-1 во вложенном таб-баре
        with dpg.plot(label="Phase envelope (P-T)", height=520, width=-1, parent=parent):
            dpg.add_plot_legend()
            dpg.add_plot_axis(dpg.mvXAxis, label="Temperature, C")
            yax = dpg.add_plot_axis(dpg.mvYAxis, label="Pressure, bar")
            if result.get("method") == "grid":
                if sta.get("T"):
                    dpg.add_scatter_series(sta["T"], sta["P"],
                                           label="Single-phase (stable)", parent=yax)
                if uns.get("T"):
                    dpg.add_scatter_series(uns["T"], uns["P"],
                                           label="Two-phase (unstable)", parent=yax)
            else:
                if env.get("T"):
                    dpg.add_line_series(env["T"], env["P"],
                                        label="Phase envelope", parent=yax)
            self._envelope_reservoir_overlay(res, p_hi, yax)

    def _envelope_reservoir_overlay(self, res, p_hi, yax) -> None:
        """Вертикаль пластовой T + маркер Psat (если посчитан) на графике огибающей."""
        if not res or res.get("T_c") is None:
            return
        tr = res["T_c"]
        dpg.add_line_series([tr, tr], [0.0, p_hi],
                            label=f"Reservoir T={self._g(tr)}C", parent=yax)
        if res.get("P_sat") is not None:
            dpg.add_scatter_series(
                [tr], [res["P_sat"]],
                label=f"Reservoir Psat={self._g(res['P_sat'])} bar", parent=yax)

    def _render_envelope_table(self, result, parent) -> None:
        tbl = result.get("table", {})
        cols = tbl.get("columns", [])
        rows = tbl.get("rows", [])
        if not cols:
            dpg.add_text("No data.", parent=parent)
            return
        formatted_rows = [[self._fmt(value) for value in row] for row in rows]
        dpg.add_button(label="Copy table", parent=parent,
                       callback=lambda: self._copy_table(
                           cols, formatted_rows, "Phase envelope data"))
        render_readonly_table(parent, cols, formatted_rows)

    # --- всплывающее окно параметров огибающей ----------------------------

    def _on_envelope_edit(self, sender, app_data, user_data) -> None:
        """«Edit parameters & run…» на вкладке — открыть диалог для узла."""
        self._open_envelope_dialog(user_data)

    def _open_envelope_dialog(self, node_id) -> None:
        """
        Модальное окно параметров огибающей перед расчётом. `node_id=None` —
        параметры нового узла (создаётся по «Run»); иначе — правка/пересчёт
        существующего узла.
        """
        composition = self._state.active_composition
        if composition is None:
            return
        node = self._state.node_by_id(node_id) if node_id is not None else None
        if node is not None:
            params = dict(node.params)
            title = "Phase envelope - edit & run"
        else:
            node_id = None
            params = pe_svc.default_envelope_params(composition)
            title = "New phase envelope"

        method = params.get("method", "ssm")
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        with dpg.window(label=title, modal=True, no_resize=True, no_collapse=True,
                        width=380, height=340,
                        pos=(max(0, w // 2 - 190), max(0, h // 2 - 170))) as win:
            self._track_modal(win)
            self._env_dialog_win = win
            self._env_dialog_node = node_id
            ids: dict = {}
            ids["method"] = dpg.add_combo(
                items=[_ENV_METHOD_SSM, _ENV_METHOD_GRID], width=220, label="Method",
                default_value=(_ENV_METHOD_GRID if method == "grid" else _ENV_METHOD_SSM),
                callback=self._on_env_method_changed)
            dpg.add_separator()

            # SSM: марш по температуре (кривая)
            with dpg.group(show=(method == "ssm")) as ssm_grp:
                dpg.add_text("Temperature march + pressure cap:")
                ids["t_min_c"] = dpg.add_input_float(
                    label="T min, C", width=150, step=0,
                    default_value=float(params.get("t_min_c", -100.0)))
                ids["t_max_c"] = dpg.add_input_float(
                    label="T max, C", width=150, step=0,
                    default_value=float(params.get("t_max_c", 200.0)))
                ids["t_step_c"] = dpg.add_input_float(
                    label="T step, C", width=150, step=0,
                    default_value=float(params.get("t_step_c", 10.0)))
                ids["p_max_bar"] = dpg.add_input_float(
                    label="P max, bar", width=150, step=0,
                    default_value=float(params.get("p_max_bar", 700.0)))
            ids["_ssm_grp"] = ssm_grp

            # Grid: скан стабильности P×T (T идёт от 0 °C)
            with dpg.group(show=(method == "grid")) as grid_grp:
                dpg.add_text("Stability scan grid (T from 0 C, P from 1 bar):")
                ids["grid_t_max_c"] = dpg.add_input_float(
                    label="T max, C", width=150, step=0,
                    default_value=float(params.get("grid_t_max_c", 300.0)))
                ids["grid_p_max_bar"] = dpg.add_input_float(
                    label="P max, bar", width=150, step=0,
                    default_value=float(params.get("grid_p_max_bar", 700.0)))
                ids["grid_t_points"] = dpg.add_input_int(
                    label="T points", width=150, step=0,
                    default_value=int(params.get("grid_t_points", 30)))
                ids["grid_p_points"] = dpg.add_input_int(
                    label="P points", width=150, step=0,
                    default_value=int(params.get("grid_p_points", 30)))
            ids["_grid_grp"] = grid_grp

            self._env_dialog_ids = ids
            dpg.add_separator()
            tres = params.get("T_res_c")
            if tres is not None:
                dpg.add_text(f"Reservoir T = {self._g(tres)} C (Psat marker)")
            dpg.add_spacer(height=4)
            with dpg.group(horizontal=True):
                dpg.add_button(label="Run", width=110,
                               callback=self._on_envelope_dialog_run)
                dpg.add_button(label="Cancel", width=110,
                               callback=lambda: dpg.delete_item(win))

    def _on_env_method_changed(self, sender, app_data) -> None:
        """Переключает видимые поля диалога под выбранный метод."""
        method = _env_method_from_label(app_data)
        ids = self._env_dialog_ids
        if dpg.does_item_exist(ids.get("_ssm_grp", 0)):
            dpg.configure_item(ids["_ssm_grp"], show=(method == "ssm"))
        if dpg.does_item_exist(ids.get("_grid_grp", 0)):
            dpg.configure_item(ids["_grid_grp"], show=(method == "grid"))

    def _on_envelope_dialog_run(self, sender, app_data, user_data) -> None:
        if self._jobs.busy:
            self._set_status("Another calculation is already running.")
            return
        composition = self._state.active_composition
        if composition is None:
            return
        ids = self._env_dialog_ids
        node_id = self._env_dialog_node
        method = _env_method_from_label(dpg.get_value(ids["method"]))
        # стартуем от полного набора параметров (сохраняем поля второго метода)
        node = self._state.node_by_id(node_id) if node_id is not None else None
        if node_id is not None and node is None:
            node_id = None
        base = (pe_svc.default_envelope_params(composition) if node is None
                else dict(node.params))
        params = dict(base)
        params["method"] = method
        try:
            if method == "grid":
                params["grid_t_max_c"] = float(dpg.get_value(ids["grid_t_max_c"]))
                params["grid_p_max_bar"] = float(dpg.get_value(ids["grid_p_max_bar"]))
                params["grid_t_points"] = int(dpg.get_value(ids["grid_t_points"]))
                params["grid_p_points"] = int(dpg.get_value(ids["grid_p_points"]))
                if params["grid_t_max_c"] <= 0 or params["grid_p_max_bar"] <= 0:
                    raise ValueError("T max and P max must be positive")
                if params["grid_t_points"] < 2 or params["grid_p_points"] < 2:
                    raise ValueError("need at least 2 grid points per axis")
            else:
                params["t_min_c"] = float(dpg.get_value(ids["t_min_c"]))
                params["t_max_c"] = float(dpg.get_value(ids["t_max_c"]))
                params["t_step_c"] = float(dpg.get_value(ids["t_step_c"]))
                params["p_max_bar"] = float(dpg.get_value(ids["p_max_bar"]))
                if params["t_max_c"] <= params["t_min_c"]:
                    raise ValueError("T max must be greater than T min")
                if params["t_step_c"] <= 0:
                    raise ValueError("T step must be positive")
                if params["p_max_bar"] <= 0:
                    raise ValueError("P max must be positive")
        except (ValueError, KeyError) as exc:
            self._set_status(f"Invalid input: {exc}")
            return

        if node_id is None:
            node_id = self._state.new_envelope(params)
        else:
            self._state.update_node_params(node_id, params, notify=False)
        assert node_id is not None

        if dpg.does_item_exist(self._env_dialog_win):
            dpg.delete_item(self._env_dialog_win)
        self._start_envelope_job(node_id)

    def _start_envelope_job(self, nid: str) -> None:
        """Запускает расчёт огибающей узла `nid` в воркер-потоке (прогресс/отмена)."""
        node = self._state.node_by_id(nid)
        composition = self._state.active_composition
        if node is None or composition is None:
            return
        params = dict(node.params)

        ref = self._state.node_ref(nid)
        if ref is None:
            return
        self._jobs.start(
            ref, "phase envelope",
            lambda: pe_svc.run_envelope(composition, params),
        )
        self._state.set_node_running(ref)
        self._set_status("Phase envelope running in parallel (cancellation unavailable)...")
        self._arm_flash_poll()
