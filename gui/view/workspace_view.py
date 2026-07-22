"""IDE-дерево модели и диспетчер вкладок рабочего пространства."""

import dearpygui.dearpygui as dpg

from gui.app_state import (
    GraphNode,
    NodeKind,
    NodeRef,
    NodeStatus,
    StateChange,
    StateChangeKind,
)
from gui.services import experiment_service as exp_svc
from gui.view.contracts import ContextBoundView

_MODEL_TREE = "model_tree"
_WORKSPACE = "workspace_content"


class WorkspaceViewMixin(ContextBoundView):
    """Навигационное дерево, context menus и маршрутизация вкладок."""

    _expanded_models: set[str]
    _expanded_cats: set[str]
    _compare_selection: list[NodeRef]
    _tabbar_id: int | None
    _tab_ids: dict[str, int]
    _tab_content_ids: dict[str, int]
    _workspace_crumb_id: int | None
    _rendered_workspace_ref: tuple[str, str] | None
    _bip_ids: dict[tuple[int, int], int]
    _flash_input_ids: dict[str, tuple[int, int]]
    _exp_input_ids: dict[str, dict[str, int]]
    _exp_chart_holder: dict[str, int]
    _lab_data_holder: dict[str, int]
    _lab_data_controls: dict[str, tuple[int, int, int]]

    def _compare_ref(self, model_id: str, node_id: str) -> NodeRef:
        model = self._state.models.get(model_id)
        variant_id = next(iter(model.variants), "base") if model else "base"
        return NodeRef(model_id, variant_id, node_id)

    def _selected_compare_nodes(self) -> list[tuple[NodeRef, GraphNode]]:
        """Живые выбранные узлы с полным межмодельным адресом."""
        out = []
        for ref in self._compare_selection:
            node = self._state.node_by_ref(ref)
            if node is not None:
                out.append((ref, node))
        return out

    def _compare_selection_count(self, kind: NodeKind,
                                 experiment_kind: str | None = None) -> int:
        return len([
            node for _ref, node in self._selected_compare_nodes()
            if node.kind is kind and (
                experiment_kind is None
                or node.params.get("kind") == experiment_kind
            )
        ])

    def _render_tree(self) -> None:
        """Показывает модели текущего проекта и узлы активной модели."""
        if not dpg.does_item_exist(_MODEL_TREE):
            return
        dpg.delete_item(_MODEL_TREE, children_only=True)

        project_id = self._state.active_project_id
        models = [model for model in self._state.models.values()
                  if project_id is None or model.project_id == project_id]
        if not models:
            dpg.add_text("No models in this project.", parent=_MODEL_TREE)
            return
        project = self._state.projects.get(project_id) if project_id else None
        project_label = project.title if project else (project_id or "Current")
        dpg.add_text(f"Models in project: {project_label}", parent=_MODEL_TREE)
        for model in models:
            expanded = model.model_id in self._expanded_models
            arrow = "v " if expanded else "> "
            # ASCII-маркер не зависит от шрифта DPG: Unicode `●` в некоторых
            # окружениях отображался как знак вопроса.
            active_mark = "[A] " if model.model_id == self._state.active_model_id else "    "
            sel = dpg.add_selectable(
                label=(f"{active_mark}{arrow}{model.title}{'  *' if model.dirty else ''}  "
                       f"[{model.n_components}c, {model.eos}]"),
                parent=_MODEL_TREE,
                default_value=(model.model_id == self._state.active_model_id),
                user_data=model.model_id, callback=self._on_model_row,
            )
            self._attach_model_context_menu(sel, model.model_id)
            # Раскрываем узлы каждой раскрытой и уже загруженной модели. Активна
            # при этом только одна модель, но несколько корней могут оставаться
            # раскрытыми одновременно — как в обычном IDE-дереве.
            if expanded and model.loaded:
                self._render_model_children(model)

    def _render_model_children(self, model) -> None:
        variant = model.variants.get("base")
        if variant is None:
            return
        is_active_model = model.model_id == self._state.active_model_id
        active_nid = variant.active_node_id if is_active_model else None

        # Composition
        comp = variant.nodes.get("composition")
        stale = "  *" if (comp and comp.status is NodeStatus.STALE) else ""
        dpg.add_selectable(
            label=f"    Composition{stale}", parent=_MODEL_TREE,
            default_value=(active_nid == "composition"),
            user_data=(model.model_id, "composition"),
            callback=self._on_tree_open_node,
        )

        # Flash (категория с историей запусков)
        cat_key = f"{model.model_id}:flash"
        cat_exp = cat_key in self._expanded_cats
        runs = variant.flash_runs()
        dpg.add_selectable(
            label=f"  {'v' if cat_exp else '>'} Flash ({len(runs)})",
            parent=_MODEL_TREE, user_data=cat_key, callback=self._on_cat_toggle,
        )
        if cat_exp:
            for run in runs:
                mark = ("[*] " if self._compare_ref(model.model_id, run.node_id)
                        in self._compare_selection else "")
                sel = dpg.add_selectable(
                    label="      " + mark + self._flash_tree_label(run),
                    parent=_MODEL_TREE, default_value=(active_nid == run.node_id),
                    user_data=(model.model_id, run.node_id),
                    callback=self._on_tree_open_node,
                )
                self._attach_flash_context_menu(sel, model.model_id, run.node_id)
            dpg.add_selectable(label="      + New flash", parent=_MODEL_TREE,
                               user_data=model.model_id, callback=self._on_new_flash)
            n_sel = self._compare_selection_count(NodeKind.FLASH)
            if n_sel >= 2:
                dpg.add_selectable(label=f"      = Compare selected ({n_sel})",
                                   parent=_MODEL_TREE, user_data=model.model_id,
                                   callback=self._on_open_compare)

        # Experiments — вложенность по типам: Experiments → DLE/CCE/Separator → запуски
        mid = model.model_id
        ecat_key = f"{mid}:exp"
        ecat_exp = ecat_key in self._expanded_cats
        exp_runs = variant.experiment_runs()
        dpg.add_selectable(
            label=f"  {'v' if ecat_exp else '>'} Experiments ({len(exp_runs)})",
            parent=_MODEL_TREE, user_data=ecat_key, callback=self._on_cat_toggle)
        if ecat_exp:
            for kind in ("cce", "dle", "separator"):
                lbl = exp_svc.EXPERIMENT_TYPES[kind]["label"]
                kruns = [r for r in exp_runs if r.params.get("kind") == kind]
                kkey = f"{mid}:exp:{kind}"
                kexp = kkey in self._expanded_cats
                dpg.add_selectable(
                    label=f"    {'v' if kexp else '>'} {lbl} ({len(kruns)})",
                    parent=_MODEL_TREE, user_data=kkey, callback=self._on_cat_toggle)
                if kexp:
                    for run in kruns:
                        mark = ("[*] " if self._compare_ref(mid, run.node_id)
                                in self._compare_selection else "")
                        sel = dpg.add_selectable(
                            label="        " + mark + self._exp_run_leaf_label(run),
                            parent=_MODEL_TREE,
                            default_value=(active_nid == run.node_id),
                            user_data=(mid, run.node_id),
                            callback=self._on_tree_open_node)
                        self._attach_exp_context_menu(sel, mid, run.node_id)
                    selected = self._compare_selection_count(
                        NodeKind.EXPERIMENT, kind,
                    )
                    if selected >= 2:
                        dpg.add_selectable(
                            label=f"        = Compare selected {lbl} ({selected})",
                            parent=_MODEL_TREE, user_data=mid,
                            callback=self._on_open_compare,
                        )
                    dpg.add_selectable(label=f"        + New {lbl}",
                                       parent=_MODEL_TREE, user_data=(mid, kind),
                                       callback=self._on_new_experiment)

        # Phase envelope — категория с историей запусков (P-T огибающая + Psat)
        pcat_key = f"{mid}:env"
        pcat_exp = pcat_key in self._expanded_cats
        env_runs = variant.envelope_runs()
        dpg.add_selectable(
            label=f"  {'v' if pcat_exp else '>'} Phase envelope ({len(env_runs)})",
            parent=_MODEL_TREE, user_data=pcat_key, callback=self._on_cat_toggle)
        if pcat_exp:
            for run in env_runs:
                sel = dpg.add_selectable(
                    label="      " + self._env_run_leaf_label(run),
                    parent=_MODEL_TREE,
                    default_value=(active_nid == run.node_id),
                    user_data=(mid, run.node_id),
                    callback=self._on_tree_open_node)
                self._attach_env_context_menu(sel, mid, run.node_id)
            dpg.add_selectable(label="      + New envelope", parent=_MODEL_TREE,
                               user_data=mid, callback=self._on_new_envelope)

    # --- обработчики дерева ----------------------------------------------

    def _attach_model_context_menu(self, item_id, model_id: str) -> None:
        """Контекстные действия корневого узла модели в рабочем дереве."""
        model = self._state.models.get(model_id)
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text(model.title if model else model_id)
            dpg.add_separator()
            dpg.add_menu_item(label="View composition (read-only)",
                              user_data=model_id,
                              callback=self._on_view_composition)
            dpg.add_menu_item(label="Save model", user_data=model_id,
                              callback=self._on_save_tree_model,
                              enabled=bool(model and model.loaded and model.dirty))
            dpg.add_menu_item(label="Duplicate model...", user_data=model_id,
                              callback=self._on_duplicate_model_confirm)

    def _on_save_tree_model(self, sender, app_data, user_data) -> None:
        """Сохраняет именно модель, из контекстного меню которой вызвано действие."""
        model_id = str(user_data)
        model = self._state.models.get(model_id)
        if model is None or not model.loaded:
            self._set_status("No loaded model to save.")
            return
        try:
            self._state.save_model(model_id)
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Save failed: {exc}")
            return
        self._set_status(f"Model '{model_id}' saved.")

    def _on_model_row(self, sender, app_data, user_data) -> None:
        mid = str(user_data)
        model = self._state.models.get(mid)
        if model is None:
            self._set_status(f"Model '{mid}' is no longer available.")
            return
        # Один клик по корню всегда переключает его состояние раскрытия. Ранее
        # для неактивной модели первый клик только активировал её, а схлопывание
        # происходило лишь вторым кликом — это выглядело как запаздывающие
        # стрелки.
        if mid in self._expanded_models:
            self._expanded_models.discard(mid)
        else:
            self._expanded_models.add(mid)
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid, notify=False)
            self._restore_workspace(mid)
            self._state.notify(StateChange(StateChangeKind.NAVIGATION))
            self._set_status(f"Active model: '{model.title}'.")
            return
        self._render_tree()

    def _on_cat_toggle(self, sender, app_data, user_data) -> None:
        cat_key = str(user_data)
        model_id, separator, _category = cat_key.partition(":")
        if separator and model_id in self._state.models:
            if model_id != self._state.active_model_id:
                self._state.set_active_model(model_id, notify=False)
                self._restore_workspace(model_id)
                self._state.notify(StateChange(StateChangeKind.NAVIGATION))
                model = self._state.models[model_id]
                self._set_status(f"Active model: '{model.title}'.")
        if cat_key in self._expanded_cats:
            self._expanded_cats.discard(cat_key)
        else:
            self._expanded_cats.add(cat_key)
        self._render_tree()

    def _on_tree_open_node(self, sender, app_data, user_data) -> None:
        mid, nid = user_data
        ref = self._compare_ref(mid, nid)
        node = self._state.node_by_ref(ref)
        # Ctrl+клик по расчёту — переключить участие в сравнении
        if (self._ctrl_down() and node is not None
                and node.kind in (NodeKind.FLASH, NodeKind.EXPERIMENT)):
            self._toggle_compare(mid, nid)
            return
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid, notify=False)
            self._restore_workspace(mid)
        self._state.open_node(nid)

    def _toggle_compare(self, model_id: str, node_id: str) -> None:
        ref = self._compare_ref(model_id, node_id)
        if ref in self._compare_selection:
            self._compare_selection.remove(ref)
        else:
            node = self._state.node_by_ref(ref)
            if node is None:
                return
            def compatible(selected_ref: NodeRef) -> bool:
                selected_node = self._state.node_by_ref(selected_ref)
                if selected_node is None or selected_node.kind is not node.kind:
                    return False
                return (node.kind is NodeKind.FLASH
                        or selected_node.params.get("kind") == node.params.get("kind"))
            self._compare_selection = [
                selected for selected in self._compare_selection if compatible(selected)
            ]
            self._compare_selection.append(ref)
        self._render_tree()
        model_count = len({selected.model_id for selected in self._compare_selection})
        self._set_status(
            f"Compare selection: {len(self._compare_selection)} run(s) "
            f"from {model_count} model(s).",
        )

    def _on_toggle_compare(self, sender, app_data, user_data) -> None:
        if isinstance(user_data, tuple) and len(user_data) == 2:
            model_id, node_id = user_data
        else:  # обратная совместимость старых callback-данных
            model_id, node_id = self._state.active_model_id, user_data
        if model_id is not None:
            self._toggle_compare(str(model_id), str(node_id))

    def _on_open_compare(self, sender, app_data, user_data) -> None:
        mid = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid, notify=False)
            self._restore_workspace(mid)
        self._state.open_compare(list(self._compare_selection))

    def _on_new_flash(self, sender, app_data, user_data) -> None:
        mid = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid, notify=False)
            self._restore_workspace(mid)
        self._state.new_flash_run()
        self._set_status("New flash tab opened - set P/T and Run.")

    def _attach_flash_context_menu(self, item_id, model_id: str,
                                   node_id: str) -> None:
        """Правый клик по листу флэша: переименовать / дублировать / удалить."""
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text("Flash run")
            dpg.add_separator()
            dpg.add_input_text(hint="rename + Enter", width=200, on_enter=True,
                               user_data=(model_id, node_id),
                               callback=self._on_flash_rename)
            dpg.add_button(label="Add / remove from compare", width=200,
                           user_data=(model_id, node_id),
                           callback=self._on_toggle_compare)
            dpg.add_button(label="Duplicate", width=200,
                           user_data=(model_id, node_id),
                           callback=self._on_flash_duplicate)
            dpg.add_button(label="Delete", width=200,
                           user_data=(model_id, node_id),
                           callback=self._on_flash_delete)

    def _node_target(self, user_data) -> tuple[str | None, str]:
        """Разбирает callback-адрес узла с обратной совместимостью."""
        if isinstance(user_data, tuple) and len(user_data) == 2:
            return str(user_data[0]), str(user_data[1])
        return self._state.active_model_id, str(user_data)

    def _activate_node_model(self, model_id: str | None) -> None:
        if model_id and model_id != self._state.active_model_id:
            self._state.set_active_model(model_id, notify=False)
            self._restore_workspace(model_id)

    def _on_flash_rename(self, sender, app_data, user_data) -> None:
        model_id, node_id = self._node_target(user_data)
        self._activate_node_model(model_id)
        self._state.rename_node(node_id, app_data)
        self._set_status(f"Flash renamed to '{app_data}'." if app_data.strip()
                         else "Flash label cleared.")

    def _on_flash_duplicate(self, sender, app_data, user_data) -> None:
        model_id, node_id = self._node_target(user_data)
        self._activate_node_model(model_id)
        self._state.duplicate_flash(node_id)
        self._set_status("Flash run duplicated.")

    def _on_flash_delete(self, sender, app_data, user_data) -> None:
        model_id, node_id = self._node_target(user_data)
        self._activate_node_model(model_id)
        self._state.delete_node(node_id)
        self._set_status("Flash run deleted.")

    # --- эксперименты: дерево ---------------------------------------------

    def _exp_status_suffix(self, node) -> str:
        if node.status is NodeStatus.RUNNING:
            return " - running"
        if node.result is not None:
            return "" if node.status is NodeStatus.FRESH else " (stale)"
        return " - (not run)"

    def _exp_tree_label(self, node) -> str:
        """Подпись эксперимента с типом (для хлебных крошек)."""
        kind = node.params.get("kind", "")
        base = exp_svc.EXPERIMENT_TYPES.get(kind, {}).get("label", kind.upper())
        core = base + self._exp_status_suffix(node)
        label = node.params.get("label")
        return f"{label}  ({core})" if label else core

    def _exp_run_leaf_label(self, node) -> str:
        """Подпись листа под своей категорией — без повтора типа, с параметрами."""
        p = node.params
        seq = node.node_id.split("_")[-1]
        desc = f"T={self._g(p.get('T_c'))}C"
        if p.get("kind") in ("dle", "separator"):
            desc += f", Pres={self._g(p.get('P_res'))}"
        core = f"#{seq}  {desc}{self._exp_status_suffix(node)}"
        label = p.get("label")
        return f"{label}  ({core})" if label else core

    def _on_new_experiment(self, sender, app_data, user_data) -> None:
        mid, kind = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid, notify=False)
            self._restore_workspace(mid)
        composition = self._state.active_composition
        if composition is None:
            return
        pressures = exp_svc.default_pressures(composition)
        t_c = round(composition.T - 273.15, 2)
        defaults = {
            "pressures": pressures,
            "T_c": t_c,
            "P_res": max(pressures),
            "stage_temps_c": [15.0] * len(pressures),
        }
        self._state.new_experiment(kind, defaults)
        self._expanded_cats.add(f"{mid}:exp")
        self._expanded_cats.add(f"{mid}:exp:{kind}")
        self._set_status(f"New {kind.upper()} tab - set stages and Run.")

    def _attach_exp_context_menu(self, item_id, model_id: str,
                                 node_id: str) -> None:
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text("Experiment")
            dpg.add_separator()
            dpg.add_input_text(hint="rename + Enter", width=200, on_enter=True,
                               user_data=(model_id, node_id),
                               callback=self._on_flash_rename)
            dpg.add_button(label="Add / remove from compare", width=200,
                           user_data=(model_id, node_id),
                           callback=self._on_toggle_compare)
            dpg.add_button(label="Duplicate", width=200,
                           user_data=(model_id, node_id),
                           callback=self._on_exp_duplicate)
            dpg.add_button(label="Delete", width=200,
                           user_data=(model_id, node_id),
                           callback=self._on_flash_delete)

    def _on_exp_duplicate(self, sender, app_data, user_data) -> None:
        model_id, node_id = self._node_target(user_data)
        self._activate_node_model(model_id)
        node = self._state.node_by_id(node_id)
        if node is not None:
            defaults = {k: v for k, v in node.params.items()
                        if k not in ("kind", "label")}
            self._state.new_experiment(node.params["kind"], defaults)
            self._set_status("Experiment duplicated.")

    # --- фазовая огибающая: дерево ----------------------------------------

    def _env_run_leaf_label(self, node) -> str:
        """Подпись листа огибающей — метод, диапазон T и статус."""
        p = node.params
        seq = node.node_id.split("_")[-1]
        if p.get("method") == "grid":
            desc = f"grid T=[0..{self._g(p.get('grid_t_max_c'))}]C"
        else:
            desc = f"ssm T=[{self._g(p.get('t_min_c'))}..{self._g(p.get('t_max_c'))}]C"
        core = f"#{seq}  {desc}{self._exp_status_suffix(node)}"
        label = p.get("label")
        return f"{label}  ({core})" if label else core

    def _on_new_envelope(self, sender, app_data, user_data) -> None:
        mid = user_data
        switched = mid != self._state.active_model_id
        if switched:
            self._state.set_active_model(mid, notify=False)
            self._restore_workspace(mid)
        if self._state.active_composition is None:
            return
        self._expanded_cats.add(f"{mid}:env")
        if switched:
            self._state.notify(StateChange(StateChangeKind.NAVIGATION))
        # параметры собираем во всплывающем окне; узел заводится по «Run»
        self._open_envelope_dialog(None)

    def _attach_env_context_menu(self, item_id, model_id: str,
                                 node_id: str) -> None:
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text("Phase envelope")
            dpg.add_separator()
            dpg.add_input_text(hint="rename + Enter", width=200, on_enter=True,
                               user_data=(model_id, node_id),
                               callback=self._on_flash_rename)
            dpg.add_button(label="Duplicate", width=200,
                           user_data=(model_id, node_id),
                           callback=self._on_env_duplicate)
            dpg.add_button(label="Delete", width=200,
                           user_data=(model_id, node_id),
                           callback=self._on_flash_delete)

    def _on_env_duplicate(self, sender, app_data, user_data) -> None:
        model_id, node_id = self._node_target(user_data)
        self._activate_node_model(model_id)
        node = self._state.node_by_id(node_id)
        if node is not None:
            defaults = {k: v for k, v in node.params.items() if k != "label"}
            self._state.new_envelope(defaults)
            self._set_status("Phase envelope duplicated.")

    # ==================================================================
    #  Рабочая область (вкладки)
    # ==================================================================

    def _render_workspace(self) -> None:
        if not dpg.does_item_exist(_WORKSPACE):
            return
        current_workspace_ref = (
            (self._state.active_model_id, self._state.active_variant_id)
            if (self._state.active_model_id is not None
                and self._state.active_variant_id is not None)
            else None
        )
        # синхронизировать закрытые крестиком вкладки до пересборки (иначе
        # они «воскресли» бы из состояния при следующей отрисовке). Tab ids
        # предыдущей модели нельзя применять к новой: exp_1 там другой узел.
        if current_workspace_ref == self._rendered_workspace_ref:
            self._reconcile_closed_tabs()
        dpg.delete_item(_WORKSPACE, children_only=True)
        self._rendered_workspace_ref = current_workspace_ref
        self._flash_input_ids = {}
        self._exp_input_ids = {}
        self._exp_chart_holder = {}
        self._lab_data_holder = {}
        self._lab_data_controls = {}

        model = self._state.active_model
        variant = self._state.active_variant
        if model is None or variant is None:
            dpg.add_text("Select a model in the tree on the left.", parent=_WORKSPACE)
            return

        # хлебные крошки
        node = self._state.active_node
        crumb = model.title + (f"   >   {self._node_crumb(node)}" if node else "")
        self._workspace_crumb_id = dpg.add_text(crumb, parent=_WORKSPACE)
        dpg.add_separator(parent=_WORKSPACE)

        if not variant.open_node_ids:
            dpg.add_text("Open a node from the tree on the left.", parent=_WORKSPACE)
            return

        self._tab_ids = {}
        self._tab_content_ids = {}
        with dpg.tab_bar(parent=_WORKSPACE, reorderable=True,
                         callback=self._on_tab_changed) as tabbar:
            self._tabbar_id = tabbar
            for nid in variant.open_node_ids:
                n = variant.nodes.get(nid)
                if n is None:
                    continue
                # closable=True -> крестик «×» рядом с названием вкладки
                with dpg.tab(label=self._tab_label(n), closable=True) as tab_id:
                    self._tab_ids[nid] = tab_id
                    page = dpg.add_group()
                    self._tab_content_ids[nid] = page
                    dpg.add_spacer(height=2, parent=page)
                    if n.kind is NodeKind.COMPOSITION:
                        self._render_composition_tab(page, n)
                    elif n.kind is NodeKind.FLASH:
                        self._render_flash_tab(page, n)
                    elif n.kind is NodeKind.EXPERIMENT:
                        self._render_experiment_tab(page, n)
                    elif n.kind is NodeKind.PHASE_ENVELOPE:
                        self._render_envelope_tab(page, n)
                    elif n.kind is NodeKind.COMPARE:
                        self._render_compare_tab(page, n)

        # синхронизировать активную вкладку
        active = variant.active_node_id
        if active in self._tab_ids:
            dpg.set_value(self._tabbar_id, self._tab_ids[active])

    def _render_node_content(self, node_id: str) -> bool:
        """Пересобирает только body одной уже открытой вкладки."""
        variant = self._state.active_variant
        node = variant.nodes.get(node_id) if variant is not None else None
        page = self._tab_content_ids.get(node_id)
        if node is None or page is None or not dpg.does_item_exist(page):
            return False

        dpg.delete_item(page, children_only=True)
        dpg.add_spacer(height=2, parent=page)
        if node.kind is NodeKind.COMPOSITION:
            self._bip_ids = {}
        elif node.kind is NodeKind.FLASH:
            self._flash_input_ids.pop(node_id, None)
        elif node.kind is NodeKind.EXPERIMENT:
            self._exp_input_ids.pop(node_id, None)
            self._exp_chart_holder.pop(node_id, None)
            self._lab_data_holder.pop(node_id, None)
            self._lab_data_controls.pop(node_id, None)

        if node.kind is NodeKind.COMPOSITION:
            self._render_composition_tab(page, node)
        elif node.kind is NodeKind.FLASH:
            self._render_flash_tab(page, node)
        elif node.kind is NodeKind.EXPERIMENT:
            self._render_experiment_tab(page, node)
        elif node.kind is NodeKind.PHASE_ENVELOPE:
            self._render_envelope_tab(page, node)
        elif node.kind is NodeKind.COMPARE:
            self._render_compare_tab(page, node)

        tab_id = self._tab_ids.get(node_id)
        if tab_id is not None and dpg.does_item_exist(tab_id):
            dpg.configure_item(tab_id, label=self._tab_label(node))
        if (self._workspace_crumb_id is not None
                and self._state.active_node is not None
                and node_id == self._state.active_node.node_id):
            model = self._state.active_model
            if model is not None:
                dpg.set_value(self._workspace_crumb_id,
                              f"{model.title}   >   {self._node_crumb(node)}")
        return True

    def _on_tab_changed(self, sender, app_data, user_data=None) -> None:
        # закрытие вкладки крестиком меняет выбор -> ловим его здесь же
        if self._reconcile_closed_tabs():
            self._render_workspace()
            return
        # app_data — id выбранной вкладки; обратное сопоставление к node_id
        nid = next((k for k, v in self._tab_ids.items() if v == app_data), None)
        if nid is not None:
            self._state.focus_node(nid)
            self._render_tree()  # обновить подсветку без пересборки вкладок

    def _reconcile_closed_tabs(self) -> bool:
        """
        Синхронизирует состояние с вкладками, закрытыми крестиком (DPG ставит
        `show=False` или удаляет элемент). Правит `open_node_ids`/`active_node_id`
        БЕЗ уведомления (вызывается из рендера/коллбэка вкладок). Возвращает
        True, если что-то закрылось.
        """
        variant = self._state.active_variant
        if variant is None:
            return False
        closed = []
        for nid, tab_id in self._tab_ids.items():
            if nid not in variant.open_node_ids:
                continue
            if (not dpg.does_item_exist(tab_id)
                    or not dpg.get_item_configuration(tab_id).get("show", True)):
                closed.append(nid)
        for nid in closed:
            self._state.close_node(nid, notify=False)
        return bool(closed)

    def _tab_label(self, node) -> str:
        if node.kind is NodeKind.COMPOSITION:
            return "Composition"
        if node.kind is NodeKind.COMPARE:
            members = node.params.get("member_refs", node.params.get("members", []))
            return f"Compare ({len(members)})"
        if node.kind is NodeKind.EXPERIMENT:
            label = node.params.get("label")
            base = exp_svc.EXPERIMENT_TYPES.get(node.params.get("kind"), {}).get(
                "label", "Exp")
            return label or f"{base} {node.node_id.split('_')[-1]}"
        if node.kind is NodeKind.PHASE_ENVELOPE:
            return node.params.get("label") or f"Envelope {node.node_id.split('_')[-1]}"
        if node.kind is NodeKind.FLASH:
            label = node.params.get("label")
            if label:
                return label
            p, t = node.params.get("P"), node.params.get("T")
            return f"Flash {self._g(p)}/{self._g(t)}"
        return node.title

    def _node_crumb(self, node) -> str:
        if node.kind is NodeKind.COMPOSITION:
            return "Composition"
        if node.kind is NodeKind.COMPARE:
            members = node.params.get("member_refs", node.params.get("members", []))
            return f"Compare ({len(members)} runs)"
        if node.kind is NodeKind.EXPERIMENT:
            return "Experiments  >  " + self._exp_tree_label(node)
        if node.kind is NodeKind.PHASE_ENVELOPE:
            return "Phase envelope  >  " + self._env_run_leaf_label(node)
        if node.kind is NodeKind.FLASH:
            return "Flash  >  " + self._flash_tree_label(node)
        return node.title

    def _flash_tree_label(self, node) -> str:
        p, t = node.params.get("P"), node.params.get("T")
        if node.status is NodeStatus.RUNNING:
            core = f"{self._g(p)}/{self._g(t)} - running"
        elif node.result is not None:
            if node.result.is_two_phase:
                phase = "2-phase"
            else:
                phase = f"1-phase {node.result.phase_type or 'unresolved'}"
            suffix = "" if node.status is NodeStatus.FRESH else " (stale)"
            core = f"{self._g(p)} bar / {self._g(t)} C - {phase}{suffix}"
        else:
            core = f"{self._g(p)} bar / {self._g(t)} C - (not run)"
        label = node.params.get("label")
        return f"{label}  ({core})" if label else core
