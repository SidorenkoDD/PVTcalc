"""IDE-дерево модели и диспетчер вкладок рабочего пространства."""

import dearpygui.dearpygui as dpg

from gui.app_state import NodeKind, NodeStatus
from gui.services import experiment_service as exp_svc
from gui.view.contracts import ContextBoundView

_MODEL_TREE = "model_tree"
_WORKSPACE = "workspace_content"


class WorkspaceViewMixin(ContextBoundView):
    """Навигационное дерево, context menus и маршрутизация вкладок."""

    _expanded_models: set[str]
    _expanded_cats: set[str]
    _compare_selection: set[str]
    _tabbar_id: int | None
    _tab_ids: dict[str, int]
    _flash_input_ids: dict[str, tuple[int, int]]
    _exp_input_ids: dict[str, dict[str, int]]
    _exp_chart_holder: dict[str, int]

    def _render_tree(self) -> None:
        """Дерево показывает ТОЛЬКО активную модель (выбор — на Projects)."""
        if not dpg.does_item_exist(_MODEL_TREE):
            return
        dpg.delete_item(_MODEL_TREE, children_only=True)

        model = self._state.active_model
        if model is None:
            dpg.add_text("No model selected.", parent=_MODEL_TREE)
            return
        expanded = model.model_id in self._expanded_models
        arrow = "v " if expanded else "> "
        dpg.add_selectable(
            label=(f"{arrow}{model.title}{'  *' if model.dirty else ''}  "
                   f"[{model.n_components}c, {model.eos}]"),
            parent=_MODEL_TREE, default_value=True,
            user_data=model.model_id, callback=self._on_model_row,
        )
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
                mark = "[*] " if run.node_id in self._compare_selection else ""
                sel = dpg.add_selectable(
                    label="      " + mark + self._flash_tree_label(run),
                    parent=_MODEL_TREE, default_value=(active_nid == run.node_id),
                    user_data=(model.model_id, run.node_id),
                    callback=self._on_tree_open_node,
                )
                self._attach_flash_context_menu(sel, run.node_id)
            dpg.add_selectable(label="      + New flash", parent=_MODEL_TREE,
                               user_data=model.model_id, callback=self._on_new_flash)
            n_sel = sum(1 for m in self._compare_selection if m in variant.nodes)
            if n_sel >= 2:
                dpg.add_selectable(label=f"      = Compare ({n_sel})",
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
                        sel = dpg.add_selectable(
                            label="        " + self._exp_run_leaf_label(run),
                            parent=_MODEL_TREE,
                            default_value=(active_nid == run.node_id),
                            user_data=(mid, run.node_id),
                            callback=self._on_tree_open_node)
                        self._attach_exp_context_menu(sel, run.node_id)
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
                self._attach_env_context_menu(sel, run.node_id)
            dpg.add_selectable(label="      + New envelope", parent=_MODEL_TREE,
                               user_data=mid, callback=self._on_new_envelope)

    # --- обработчики дерева ----------------------------------------------

    def _on_model_row(self, sender, app_data, user_data) -> None:
        # модель уже активна (единственная в дереве) — только разворот
        mid = user_data
        if mid in self._expanded_models:
            self._expanded_models.discard(mid)
        else:
            self._expanded_models.add(mid)
        self._render_tree()

    def _on_cat_toggle(self, sender, app_data, user_data) -> None:
        cat_key = user_data
        if cat_key in self._expanded_cats:
            self._expanded_cats.discard(cat_key)
        else:
            self._expanded_cats.add(cat_key)
        self._render_tree()

    def _on_tree_open_node(self, sender, app_data, user_data) -> None:
        mid, nid = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid)
        node = self._state.node_by_id(nid)
        # Ctrl+клик по флэшу — переключить участие в сравнении
        if self._ctrl_down() and node is not None and node.kind is NodeKind.FLASH:
            self._toggle_compare(nid)
            return
        self._state.open_node(nid)

    def _toggle_compare(self, node_id: str) -> None:
        if node_id in self._compare_selection:
            self._compare_selection.discard(node_id)
        else:
            self._compare_selection.add(node_id)
        self._render_tree()
        self._set_status(f"Compare selection: {len(self._compare_selection)} run(s).")

    def _on_toggle_compare(self, sender, app_data, user_data) -> None:
        self._toggle_compare(user_data)

    def _on_open_compare(self, sender, app_data, user_data) -> None:
        mid = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid)
        variant = self._state.active_variant
        if variant is None:
            return
        members = [m for m in self._compare_selection
                   if m in variant.nodes]
        self._state.open_compare(members)

    def _on_new_flash(self, sender, app_data, user_data) -> None:
        mid = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid)
        self._state.new_flash_run()
        self._set_status("New flash tab opened - set P/T and Run.")

    def _attach_flash_context_menu(self, item_id, node_id: str) -> None:
        """Правый клик по листу флэша: переименовать / дублировать / удалить."""
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text("Flash run")
            dpg.add_separator()
            dpg.add_input_text(hint="rename + Enter", width=200, on_enter=True,
                               user_data=node_id, callback=self._on_flash_rename)
            dpg.add_button(label="Add / remove from compare", width=200,
                           user_data=node_id, callback=self._on_toggle_compare)
            dpg.add_button(label="Duplicate", width=200, user_data=node_id,
                           callback=self._on_flash_duplicate)
            dpg.add_button(label="Delete", width=200, user_data=node_id,
                           callback=self._on_flash_delete)

    def _on_flash_rename(self, sender, app_data, user_data) -> None:
        self._state.rename_node(user_data, app_data)
        self._set_status(f"Flash renamed to '{app_data}'." if app_data.strip()
                         else "Flash label cleared.")

    def _on_flash_duplicate(self, sender, app_data, user_data) -> None:
        self._state.duplicate_flash(user_data)
        self._set_status("Flash run duplicated.")

    def _on_flash_delete(self, sender, app_data, user_data) -> None:
        self._state.delete_node(user_data)
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
            self._state.set_active_model(mid)
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

    def _attach_exp_context_menu(self, item_id, node_id: str) -> None:
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text("Experiment")
            dpg.add_separator()
            dpg.add_input_text(hint="rename + Enter", width=200, on_enter=True,
                               user_data=node_id, callback=self._on_flash_rename)
            dpg.add_button(label="Duplicate", width=200, user_data=node_id,
                           callback=self._on_exp_duplicate)
            dpg.add_button(label="Delete", width=200, user_data=node_id,
                           callback=self._on_flash_delete)

    def _on_exp_duplicate(self, sender, app_data, user_data) -> None:
        node = self._state.node_by_id(user_data)
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
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid)
        if self._state.active_composition is None:
            return
        self._expanded_cats.add(f"{mid}:env")
        # параметры собираем во всплывающем окне; узел заводится по «Run»
        self._open_envelope_dialog(None)

    def _attach_env_context_menu(self, item_id, node_id: str) -> None:
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text("Phase envelope")
            dpg.add_separator()
            dpg.add_input_text(hint="rename + Enter", width=200, on_enter=True,
                               user_data=node_id, callback=self._on_flash_rename)
            dpg.add_button(label="Duplicate", width=200, user_data=node_id,
                           callback=self._on_env_duplicate)
            dpg.add_button(label="Delete", width=200, user_data=node_id,
                           callback=self._on_flash_delete)

    def _on_env_duplicate(self, sender, app_data, user_data) -> None:
        node = self._state.node_by_id(user_data)
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
        # синхронизировать закрытые крестиком вкладки до пересборки (иначе
        # они «воскресли» бы из состояния при следующей отрисовке)
        self._reconcile_closed_tabs()
        dpg.delete_item(_WORKSPACE, children_only=True)
        self._flash_input_ids = {}
        self._exp_input_ids = {}
        self._exp_chart_holder = {}

        model = self._state.active_model
        variant = self._state.active_variant
        if model is None or variant is None:
            dpg.add_text("Select a model in the tree on the left.", parent=_WORKSPACE)
            return

        # хлебные крошки
        node = self._state.active_node
        crumb = model.title + (f"   >   {self._node_crumb(node)}" if node else "")
        dpg.add_text(crumb, parent=_WORKSPACE)
        dpg.add_separator(parent=_WORKSPACE)

        if not variant.open_node_ids:
            dpg.add_text("Open a node from the tree on the left.", parent=_WORKSPACE)
            return

        self._tab_ids = {}
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
            return f"Compare ({len(node.params.get('members', []))})"
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
            return f"Compare ({len(node.params.get('members', []))} runs)"
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
