"""Редактор состава: свойства компонентов, корреляции и BIP."""

import dearpygui.dearpygui as dpg

from gui.app_state import NodeStatus
from gui.services import composition_service as comp_svc
from gui.view.contracts import ContextBoundView

_BIP_CELL_W = 60

_COMPOSITION_COLUMNS: list[tuple[str, str]] = [
    ("molar_mass", "M, g/mol"),
    ("critical_temperature", "Tc, K"),
    ("critical_pressure", "Pc, bar"),
    ("acentric_factor", "omega"),
    ("critical_volume", "Vc"),
    ("shift_parameter", "shift"),
    ("Kw", "Kw"),
]


class CompositionViewMixin(ContextBoundView):
    """Рендер и callbacks вкладки Composition."""

    _bip_ids: dict[tuple[int, int], int]

    def _render_composition_tab(self, parent, node) -> None:
        composition = self._state.active_composition
        if composition is None:
            dpg.add_text("No composition.", parent=parent)
            return

        status_line = dpg.add_text(f"Node status: {node.status.name}", parent=parent)
        if node.status is NodeStatus.STALE:
            dpg.bind_item_theme(status_line, self._theme_stale())
        if node.error:
            err = dpg.add_text(f"Error: {node.error}", parent=parent)
            dpg.bind_item_theme(err, self._theme_stale())
        dpg.add_separator(parent=parent)

        with dpg.tab_bar(parent=parent):
            with dpg.tab(label="Properties & correlations"):
                p = dpg.add_group()
                self._render_controls(composition, node, p)
                dpg.add_separator(parent=p)
                self._render_composition_table(composition, p)
            with dpg.tab(label="BIP matrix"):
                p = dpg.add_group()
                self._render_bip_matrix(composition, p)

    def _render_controls(self, composition, node, parent) -> None:
        params = node.params if node else {}
        eos_value = params.get("eos", composition.eos_name.value)
        correlations = params.get("correlations", comp_svc.default_correlations())

        with dpg.group(horizontal=True, parent=parent):
            dpg.add_combo(items=comp_svc.EOS_OPTIONS, default_value=eos_value,
                          width=140, label="EOS", callback=self._on_eos_changed)
            dpg.add_button(label="Recalculate properties",
                           callback=lambda: self._state.recalculate_composition())

        if comp_svc.has_c7_plus(composition):
            with dpg.tree_node(label="C7+ correlations", parent=parent,
                               default_open=True):
                for prop, prop_label in comp_svc.CORRELATION_LABELS:
                    dpg.add_combo(items=comp_svc.CORRELATION_OPTIONS[prop],
                                  default_value=correlations.get(prop, ""), width=220,
                                  label=prop_label, user_data=prop,
                                  callback=self._on_correlation_changed)
        else:
            dpg.add_text("No C7+ components - correlations not used.", parent=parent)

    def _render_composition_table(self, composition, parent) -> None:
        comp = composition.composition
        data = composition.composition_data

        with dpg.group(horizontal=True, parent=parent):
            dpg.add_text(f"Sum zi = {comp_svc.sum_zi(composition):.6f}")
            dpg.add_button(label="Normalize zi",
                           callback=lambda: self._state.normalize_composition())

        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, height=-1, freeze_rows=1):
            dpg.add_table_column(label="Component")
            dpg.add_table_column(label="zi, frac.")
            for _key, title in _COMPOSITION_COLUMNS:
                dpg.add_table_column(label=title)
            for name in comp:
                with dpg.table_row():
                    dpg.add_text(name)
                    dpg.add_input_float(default_value=float(comp[name]), width=90,
                                        step=0, format="%.5f", on_enter=True,
                                        user_data=name, callback=self._on_zi_edited)
                    for key, _title in _COMPOSITION_COLUMNS:
                        value = data.get(key, {}).get(name)
                        if value is None:
                            dpg.add_text("-")
                        else:
                            dpg.add_input_float(default_value=float(value), width=90,
                                                step=0, format="%.5g", on_enter=True,
                                                user_data=(name, key),
                                                callback=self._on_property_edited)

    def _render_bip_matrix(self, composition, parent) -> None:
        names = list(composition.composition.keys())
        data = composition.composition_data.get("bip", {})
        self._bip_ids = {}
        dpg.add_text("Symmetric matrix - editing cell (i, j) mirrors to (j, i). "
                     "Diagonal fixed at 0. Press Enter to apply.", parent=parent)
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       scrollX=True, scrollY=True, height=-1, freeze_rows=1,
                       freeze_columns=1, policy=dpg.mvTable_SizingFixedFit):
            dpg.add_table_column(label="")
            for name_j in names:
                dpg.add_table_column(label=name_j)
            for i, name_i in enumerate(names):
                with dpg.table_row():
                    dpg.add_text(name_i)
                    for j, name_j in enumerate(names):
                        if i == j:
                            dpg.add_input_float(default_value=0.0, width=_BIP_CELL_W,
                                                step=0, format="%.3f", readonly=True,
                                                enabled=False)
                        else:
                            value = float(data.get(name_i, {}).get(name_j, 0.0))
                            iid = dpg.add_input_float(
                                default_value=value, width=_BIP_CELL_W, step=0,
                                format="%.3f", on_enter=True, user_data=(i, j),
                                callback=self._on_bip_cell_edited)
                            self._bip_ids[(i, j)] = iid

    # --- обработчики Composition ----------------------------------------

    def _on_eos_changed(self, sender, app_data, user_data) -> None:
        self._state.set_composition_eos(app_data)
        self._set_status(f"EOS set to {app_data}.")

    def _on_correlation_changed(self, sender, app_data, user_data) -> None:
        self._state.set_correlation(user_data, app_data)
        self._set_status(f"Correlation {user_data} = '{app_data}' (press Recalculate).")

    def _on_zi_edited(self, sender, app_data, user_data) -> None:
        try:
            self._state.edit_zi(user_data, app_data)
        except ValueError as exc:
            self._show_input_validation({"value": sender}, {"value": str(exc)})
            composition = self._state.active_composition
            if composition is not None:
                dpg.set_value(sender, composition.composition[user_data])
            self._set_status(f"Invalid mole fraction: {exc}")
            return
        self._show_input_validation({"value": sender}, {})
        self._set_status(f"zi[{user_data}] = {app_data:.5f} (not normalized).")

    def _on_property_edited(self, sender, app_data, user_data) -> None:
        name, key = user_data
        try:
            self._state.edit_component_property(name, key, app_data)
        except ValueError as exc:
            self._show_input_validation({"value": sender}, {"value": str(exc)})
            composition = self._state.active_composition
            if composition is not None:
                dpg.set_value(sender, composition.composition_data[key][name])
            self._set_status(f"Invalid property: {exc}")
            return
        self._show_input_validation({"value": sender}, {})
        self._set_status(f"{key}[{name}] = {app_data:.5g}.")

    def _on_bip_cell_edited(self, sender, app_data, user_data) -> None:
        i, j = user_data
        composition = self._state.active_composition
        if composition is None:
            return
        names = list(composition.composition.keys())
        try:
            self._state.edit_bip(names[i], names[j], app_data)
        except ValueError as exc:
            self._show_input_validation({"value": sender}, {"value": str(exc)})
            value = comp_svc.get_bip(composition, names[i], names[j])
            dpg.set_value(sender, value)
            self._set_status(f"Invalid BIP: {exc}")
            return
        self._show_input_validation({"value": sender}, {})
        mirror = self._bip_ids.get((j, i))
        if mirror is not None and dpg.does_item_exist(mirror):
            dpg.set_value(mirror, app_data)
        self._set_status(f"BIP[{names[i]},{names[j]}] = {app_data:.5f}.")
