"""Вкладки PVT-экспериментов, таблицы и графики результатов."""

from typing import Any

import dearpygui.dearpygui as dpg

from gui.app_state import NodeStatus
from gui.services import experiment_service as exp_svc


class ExperimentViewMixin:
    """Рендер и callbacks CCE/DLE/Separator."""

    _state: Any
    _jobs: Any
    _exp_input_ids: Any
    _exp_chart_holder: Any
    _fmt: Any
    _g: Any
    _set_status: Any
    _theme_stale: Any
    _arm_flash_poll: Any
    _on_flash_cancel: Any

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

    def _render_exp_table(self, result, parent, cols) -> None:
        allcols = result["columns"]
        idxs = [(c, allcols.index(c)) for c in cols if c in allcols]
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True, height=-1,
                       freeze_rows=1, freeze_columns=1):
            for c, _ in idxs:
                dpg.add_table_column(label=c)
            for row in result["rows"]:
                with dpg.table_row():
                    for _, i in idxs:
                        dpg.add_text(self._fmt(row[i]))

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
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True, height=-1,
                       freeze_rows=1, freeze_columns=1):
            dpg.add_table_column(label="Component")
            for st in stages:
                dpg.add_table_column(label=f"{self._g(st.get('pressure'))} bar")
            for name in names:
                with dpg.table_row():
                    dpg.add_text(name)
                    for st in stages:
                        d = st.get(phase)
                        dpg.add_text(self._fmt(d.get(name) if isinstance(d, dict)
                                               else None))

    def _render_exp_chart(self, result, parent, nid) -> None:
        node = self._state.node_by_id(nid)
        charts = list(result.get("charts", []))
        for c in (node.params.get("extra_charts", []) if node else []):
            if c not in charts:
                charts.append(c)
        with dpg.group(horizontal=True, parent=parent):
            dpg.add_combo(items=result.get("plot_all", []), width=240,
                          label="Add chart for column", user_data=nid,
                          callback=self._on_exp_add_chart)
        holder = dpg.add_group(parent=parent)
        self._exp_chart_holder[nid] = holder
        if not charts:
            dpg.add_text("No plottable columns.", parent=holder)
        for col in charts:
            self._add_one_chart(result, col, holder)

    @staticmethod
    def _exp_x_range(result) -> tuple[float, float] | None:
        """Диапазон оси X (давления) по данным эксперимента, с малым отступом."""
        cols = result.get("columns", [])
        x = result.get("x")
        if x not in cols:
            return None
        xi = cols.index(x)
        xv = [row[xi] for row in result["rows"] if row[xi] is not None]
        if not xv:
            return None
        lo, hi = min(xv), max(xv)
        margin = (hi - lo) * 0.02 if hi > lo else 1.0
        return lo - margin, hi + margin

    def _add_one_chart(self, result, col, parent) -> None:
        xs, ys = exp_svc.series_for_plot(result, col)
        if not xs:
            return
        xr = self._exp_x_range(result)
        with dpg.plot(label=f"{col} vs pressure", height=220, width=-1, parent=parent):
            dpg.add_plot_legend()
            xax = dpg.add_plot_axis(dpg.mvXAxis, label="Pressure, bar")
            yax = dpg.add_plot_axis(dpg.mvYAxis, label=col)
            dpg.add_line_series(xs, ys, label=col, parent=yax)
            if xr is not None:
                dpg.set_axis_limits(xax, xr[0], xr[1])

    def _on_exp_add_chart(self, sender, app_data, user_data) -> None:
        nid = user_data
        col = app_data
        node = self._state.node_by_id(nid)
        if not col or node is None or not isinstance(node.result, dict):
            return
        extra = node.params.setdefault("extra_charts", [])
        if col not in extra and col not in node.result.get("charts", []):
            extra.append(col)
        holder = self._exp_chart_holder.get(nid)
        if holder is not None and dpg.does_item_exist(holder):
            self._add_one_chart(node.result, col, holder)
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
