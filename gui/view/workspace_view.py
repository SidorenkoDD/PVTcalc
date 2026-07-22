"""IDE-дерево модели и диспетчер вкладок рабочего пространства."""

import math

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
from gui.services import lab_data_service as lab_svc
from gui.services import project_service as proj_svc
from gui.view.contracts import ContextBoundView

_MODEL_TREE = "model_tree"
_WORKSPACE = "workspace_content"


class WorkspaceViewMixin(ContextBoundView):
    """Навигационное дерево, context menus и маршрутизация вкладок."""

    _expanded_models: set[str]
    _expanded_cats: set[str]
    _tree_query: str
    _search_selected_model_ids: set[str]
    _selected_tree_model_id: str | None
    _compare_selection: list[NodeRef]
    _tabbar_id: int | None
    _tab_ids: dict[str, int]
    _tab_content_ids: dict[str, int]
    _overview_tab_id: int | None
    _workspace_crumb_id: int | None
    _rendered_workspace_ref: tuple[str, str] | None
    _bip_ids: dict[tuple[int, int], int]
    _flash_input_ids: dict[str, tuple[int, int]]
    _exp_input_ids: dict[str, dict[str, int]]
    _exp_chart_holder: dict[str, int]
    _lab_data_holder: dict[str, int]
    _lab_data_controls: dict[str, tuple[int, int, int]]
    _lab_catalog_editor: dict | None
    _lab_catalog_modal: int | None

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
        dpg.add_text(
            "Shortcuts: Ctrl+S save | Ctrl+Z/Y undo/redo | Del delete selected "
            "tree item | Esc close dialog.",
            parent=_MODEL_TREE, wrap=290,
        )
        with dpg.group(horizontal=True, parent=_MODEL_TREE):
            dpg.add_input_text(
                hint="Find calculations: DLE, Flash... + Enter", width=255,
                default_value=self._tree_query, callback=self._on_tree_query,
                on_enter=True,
            )
            if self._tree_query:
                dpg.add_button(label="Clear", small=True,
                               callback=self._on_tree_search_clear)
        if self._tree_query:
            self._render_search_results(models)
            return
        self._render_project_lab_data_tree(project_id)
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

    def _render_project_lab_data_tree(self, project_id: str | None) -> None:
        """Separate shared-data branch; model-local datasets stay out of it."""
        if not project_id:
            return
        try:
            datasets = lab_svc.list_datasets(self._state.db_path, project_id)
        except lab_svc.LabDataStoreError as exc:
            dpg.add_text(f"Project Lab Data unavailable: {exc}", parent=_MODEL_TREE,
                         wrap=290)
            return
        key = f"{project_id}:lab"
        expanded = key in self._expanded_cats
        dpg.add_selectable(
            label=f"{'v' if expanded else '>'} Project Lab Data ({len(datasets)})",
            parent=_MODEL_TREE, user_data=key, callback=self._on_cat_toggle,
        )
        if not expanded:
            return
        if not datasets:
            dpg.add_text("    Create a measured-data set for this project.",
                         parent=_MODEL_TREE, wrap=270)
        for dataset in datasets:
            label = (f"    {dataset['experiment_kind'].upper()} — {dataset['title']} "
                     f"({len(dataset['rows'])} point(s))")
            dpg.add_selectable(label=label, parent=_MODEL_TREE,
                               user_data=(dataset["dataset_id"], "project", None),
                               callback=self._on_project_lab_dataset)
        for kind in ("dle", "cce", "separator"):
            dpg.add_selectable(label=f"    + New project {kind.upper()} Lab Data",
                               parent=_MODEL_TREE,
                               user_data=(kind, "project", None),
                               callback=self._on_new_lab_dataset)

    def _on_project_lab_dataset(self, sender, app_data, user_data) -> None:
        dataset_id, scope, model_id = user_data
        project_id = self._state.active_project_id
        if not project_id:
            return
        dataset = lab_svc.get_dataset(self._state.db_path, project_id, dataset_id,
                                      model_id=model_id)
        if dataset is not None:
            self._open_lab_dataset_editor(dataset)

    def _render_model_lab_data_tree(self, model) -> None:
        datasets = [dataset for dataset in lab_svc.list_datasets(
            self._state.db_path, model.project_id, model_id=model.model_id)
            if dataset["scope"] == "model"]
        key = f"{model.model_id}:model_lab"
        expanded = key in self._expanded_cats
        dpg.add_selectable(
            label=f"  {'v' if expanded else '>'} Model Lab Data ({len(datasets)})",
            parent=_MODEL_TREE, user_data=key, callback=self._on_cat_toggle,
        )
        if not expanded:
            return
        for dataset in datasets:
            dpg.add_selectable(
                label=(f"      {dataset['experiment_kind'].upper()} — "
                       f"{dataset['title']} ({len(dataset['rows'])} point(s))"),
                parent=_MODEL_TREE,
                user_data=(dataset["dataset_id"], "model", model.model_id),
                callback=self._on_project_lab_dataset,
            )
        for kind in ("dle", "cce", "separator"):
            dpg.add_selectable(label=f"      + New {kind.upper()} Lab Data",
                               parent=_MODEL_TREE,
                               user_data=(kind, "model", model.model_id),
                               callback=self._on_new_lab_dataset)

    def _on_new_lab_dataset(self, sender, app_data, user_data) -> None:
        kind, scope, model_id = user_data
        self._open_lab_dataset_editor({
            "dataset_id": None,
            "title": f"{str(kind).upper()} Lab Data",
            "experiment_kind": str(kind),
            "scope": str(scope),
            "model_id": model_id,
            "columns": list(exp_svc.EXPERIMENT_TYPES[str(kind)]["lab_columns"]),
            "rows": [],
            "conditions": {},
        })

    def _open_lab_dataset_editor(self, dataset) -> None:
        """Opens the single manual editor used by project and model catalogs."""
        if self._lab_catalog_modal and dpg.does_item_exist(self._lab_catalog_modal):
            dpg.delete_item(self._lab_catalog_modal)
        data = dict(dataset)
        data["rows"] = [list(row) for row in data.get("rows", [])]
        data["conditions"] = dict(data.get("conditions") or {})
        self._lab_catalog_editor = data
        title = ("New " if data.get("dataset_id") is None else "Edit ") + "Lab Data"
        with dpg.window(label=title, modal=True, no_collapse=True,
                        width=850, height=620) as win:
            self._track_modal(win)
            self._lab_catalog_modal = win
            data["title_id"] = dpg.add_input_text(
                label="Name", width=430, default_value=data["title"], parent=win)
            dpg.add_text(
                f"{data['scope'].title()} scope · "
                f"{data['experiment_kind'].upper()} · manual input only",
                parent=win,
            )
            with dpg.group(horizontal=True, parent=win):
                data["t_id"] = dpg.add_input_text(
                    label="T, C (optional)", width=160,
                    default_value=(self._g(data["conditions"]["T_c"])
                                   if "T_c" in data["conditions"] else ""))
                if data["experiment_kind"] in ("dle", "separator"):
                    data["p_id"] = dpg.add_input_text(
                        label="P res, bar (optional)", width=180,
                        default_value=(self._g(data["conditions"]["P_res"])
                                       if "P_res" in data["conditions"] else ""))
            with dpg.group(horizontal=True, parent=win):
                dpg.add_button(label="Add point", callback=self._on_catalog_lab_add_row)
                dpg.add_button(label="Remove last", callback=self._on_catalog_lab_remove_row,
                               enabled=True)
                dpg.add_button(label="Save", callback=self._on_catalog_lab_save)
                dpg.add_button(label="Cancel", callback=self._on_catalog_lab_cancel)
            data["holder"] = dpg.add_group(parent=win)
            self._render_catalog_lab_table()

    def _render_catalog_lab_table(self) -> None:
        data = self._lab_catalog_editor
        if not data or not dpg.does_item_exist(data["holder"]):
            return
        dpg.delete_item(data["holder"], children_only=True)
        rows = data["rows"]
        columns = data["columns"]
        if not rows:
            dpg.add_text("No measured points. Add them manually one row at a time.",
                         parent=data["holder"])
            return
        with dpg.table(parent=data["holder"], header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True, height=350):
            dpg.add_table_column(label="#", width_fixed=True, width=35)
            for column in columns:
                dpg.add_table_column(label=column)
            for row_index, row in enumerate(rows):
                with dpg.table_row():
                    dpg.add_text(str(row_index + 1))
                    for column_index, value in enumerate(row):
                        dpg.add_input_text(
                            default_value="" if value is None else self._g(value),
                            width=115, user_data=(row_index, column_index),
                            callback=self._on_catalog_lab_cell,
                        )

    def _on_catalog_lab_add_row(self, sender=None, app_data=None, user_data=None) -> None:
        data = self._lab_catalog_editor
        if data is None:
            return
        data["rows"].append([None] * len(data["columns"]))
        self._render_catalog_lab_table()

    def _on_catalog_lab_remove_row(self, sender=None, app_data=None, user_data=None) -> None:
        data = self._lab_catalog_editor
        if data is None or not data["rows"]:
            return
        data["rows"].pop()
        self._render_catalog_lab_table()

    def _on_catalog_lab_cell(self, sender, app_data, user_data) -> None:
        data = self._lab_catalog_editor
        if data is None:
            return
        row, column = user_data
        text = str(app_data).strip().replace(",", ".")
        try:
            value = None if not text else float(text)
            if value is not None and not math.isfinite(value):
                raise ValueError
        except ValueError:
            self._set_status("Lab Data values must be finite numbers.")
            return
        data["rows"][row][column] = value

    def _catalog_lab_conditions(self, data) -> dict[str, float]:
        conditions: dict[str, float] = {}
        for key, field in (("T_c", "t_id"), ("P_res", "p_id")):
            item = data.get(field)
            raw = str(dpg.get_value(item)).strip() if item else ""
            if not raw:
                continue
            try:
                value = float(raw.replace(",", "."))
            except ValueError:
                raise ValueError(f"{key} must be a number") from None
            if not math.isfinite(value):
                raise ValueError(f"{key} must be finite")
            conditions[key] = value
        return conditions

    def _on_catalog_lab_save(self, sender=None, app_data=None, user_data=None) -> None:
        data = self._lab_catalog_editor
        project_id = self._state.active_project_id
        if data is None or not project_id:
            return
        try:
            title = str(dpg.get_value(data["title_id"])).strip()
            if not title:
                raise ValueError("Name is required")
            conditions = self._catalog_lab_conditions(data)
            if data.get("dataset_id"):
                saved = lab_svc.update_dataset(
                    self._state.db_path, project_id, data["dataset_id"], title=title,
                    columns=data["columns"], rows=data["rows"], conditions=conditions,
                    model_id=data.get("model_id"),
                )
            else:
                saved = lab_svc.create_dataset(
                    self._state.db_path, project_id, title=title,
                    experiment_kind=data["experiment_kind"], columns=data["columns"],
                    rows=data["rows"], conditions=conditions, scope=data["scope"],
                    model_id=data.get("model_id"),
                )
            if saved is None:
                raise ValueError("Lab Data dataset is no longer available")
        except (ValueError, lab_svc.LabDataStoreError) as exc:
            self._set_status(f"Could not save Lab Data: {exc}")
            return
        self._on_catalog_lab_cancel()
        self._render_tree()
        self._render_workspace()
        self._set_status("Lab Data saved.")

    def _on_catalog_lab_cancel(self, sender=None, app_data=None, user_data=None) -> None:
        if self._lab_catalog_modal and dpg.does_item_exist(self._lab_catalog_modal):
            dpg.delete_item(self._lab_catalog_modal)
        self._lab_catalog_modal = None
        self._lab_catalog_editor = None

    def _on_tree_query(self, sender, app_data, user_data=None) -> None:
        query = str(app_data).strip().lower()
        if query != self._tree_query:
            self._search_selected_model_ids.clear()
        self._tree_query = query
        self._render_tree()

    def _on_tree_search_clear(self, sender, app_data, user_data=None) -> None:
        self._tree_query = ""
        self._search_selected_model_ids.clear()
        self._render_tree()

    def _search_node_text(self, node) -> str:
        if node.kind is NodeKind.FLASH:
            label = self._flash_tree_label(node)
        elif node.kind is NodeKind.EXPERIMENT:
            label = self._exp_tree_label(node)
        elif node.kind is NodeKind.PHASE_ENVELOPE:
            label = self._env_run_leaf_label(node)
        else:
            label = node.title
        return " ".join(filter(None, (
            node.node_id, node.title, node.kind.name,
            str(node.params.get("kind") or ""), label,
        ))).lower()

    def _search_hits(self, model) -> list[tuple[str, str]]:
        """Возвращает расчёты модели, подходящие под запрос, без её загрузки."""
        if not self._tree_query:
            return []
        variant = model.variants.get("base")
        if model.loaded and variant is not None:
            return [
                (node.node_id, self._search_node_label(node))
                for node in variant.nodes.values()
                if (node.kind not in (NodeKind.COMPOSITION, NodeKind.COMPARE)
                    and self._tree_query in self._search_node_text(node))
            ]
        records = model.summary.results_brief if model.summary is not None else ()
        hits: list[tuple[str, str]] = []
        for record in records:
            node_id = record.get("node_id") if isinstance(record, dict) else None
            if not isinstance(node_id, str):
                continue
            text = " ".join(str(record.get(key) or "") for key in (
                "node_id", "module", "kind", "experiment_kind",
            )).lower()
            if self._tree_query in text:
                hits.append((node_id, self._saved_search_label(record)))
        return hits

    def _search_node_label(self, node) -> str:
        if node.kind is NodeKind.EXPERIMENT:
            kind = str(node.params.get("kind") or "experiment")
            return f"{kind.upper()} — {self._exp_run_leaf_label(node)}"
        if node.kind is NodeKind.FLASH:
            return "Flash — " + self._flash_tree_label(node)
        return "Phase envelope — " + self._env_run_leaf_label(node)

    @staticmethod
    def _saved_search_label(record: dict) -> str:
        kind = str(record.get("experiment_kind") or record.get("kind")
                   or record.get("module") or "Calculation")
        node_id = str(record.get("node_id") or "")
        return f"{kind.upper()} — {node_id} (saved)"

    def _render_search_results(self, models) -> None:
        matches = [(model, self._search_hits(model)) for model in models]
        matches = [(model, hits) for model, hits in matches if hits]
        matching_model_ids = {model.model_id for model, _hits in matches}
        self._search_selected_model_ids.intersection_update(matching_model_ids)
        total = sum(len(hits) for _model, hits in matches)
        dpg.add_text(
            f"Search results for '{self._tree_query}': {total} calculation(s) "
            f"in {len(matches)} model(s).",
            parent=_MODEL_TREE, wrap=290,
        )
        if not matches:
            dpg.add_text("No calculations match this search.", parent=_MODEL_TREE)
            return
        selected_count = len(self._search_selected_model_ids)
        dpg.add_text(
            "Select model groups, then compare the matching calculations.",
            parent=_MODEL_TREE, wrap=290,
        )
        dpg.add_button(
            label=f"Compare selected models ({selected_count})",
            parent=_MODEL_TREE, callback=self._on_compare_search_models,
            enabled=selected_count >= 2,
        )
        for model, hits in matches:
            with dpg.group(horizontal=True, parent=_MODEL_TREE):
                dpg.add_checkbox(
                    label="", default_value=model.model_id in self._search_selected_model_ids,
                    user_data=model.model_id, callback=self._on_search_model_selected,
                )
                dpg.add_selectable(
                    label=f"{model.title}  [{len(hits)} match(es)]",
                    user_data=model.model_id, callback=self._on_model_row,
                )
            for node_id, label in hits:
                dpg.add_selectable(
                    label="      " + label, parent=_MODEL_TREE,
                    user_data=(model.model_id, node_id),
                    callback=self._on_tree_open_node,
                )

    def _on_search_model_selected(self, sender, app_data, user_data) -> None:
        model_id = str(user_data)
        if bool(app_data):
            self._search_selected_model_ids.add(model_id)
        else:
            self._search_selected_model_ids.discard(model_id)
        self._render_tree()

    def _on_compare_search_models(self, sender, app_data, user_data=None) -> None:
        """Открывает Compare для совпадающих запусков отмеченных моделей."""
        project_id = self._state.active_project_id
        models = [model for model in self._state.models.values()
                  if (model.model_id in self._search_selected_model_ids
                      and (project_id is None or model.project_id == project_id))]
        if len(models) < 2 or not self._tree_query:
            self._set_status("Select at least two model groups from the search results.")
            return

        candidates: list[tuple[NodeRef, GraphNode]] = []
        for model in models:
            if model.model_id != self._state.active_model_id:
                self._state.set_active_model(model.model_id, notify=False)
            self._restore_workspace(model.model_id)
            variant = model.variants.get("base")
            if variant is None:
                continue
            candidates.extend(
                (NodeRef(model.model_id, variant.variant_id, node.node_id), node)
                for node in variant.nodes.values()
                if (node.kind in (NodeKind.FLASH, NodeKind.EXPERIMENT)
                    and self._tree_query in self._search_node_text(node))
            )

        groups: dict[tuple[NodeKind, str | None], list[NodeRef]] = {}
        for ref, node in candidates:
            key = (node.kind, str(node.params.get("kind"))
                   if node.kind is NodeKind.EXPERIMENT else None)
            groups.setdefault(key, []).append(ref)
        compatible_groups = [refs for refs in groups.values() if len(refs) >= 2]
        if len(compatible_groups) != 1:
            self._set_status(
                "Search a single calculation type (for example DLE) before Compare.",
            )
            return

        refs = compatible_groups[0]
        self._compare_selection = refs
        compare_members: list[str | NodeRef] = []
        compare_members.extend(refs)
        if self._state.open_compare(compare_members) is None:
            self._set_status("Selected calculations could not be compared.")
            return
        self._set_status(
            f"Compare opened: {len(refs)} matching calculation(s) from "
            f"{len({ref.model_id for ref in refs})} model(s).",
        )

    def _render_model_children(self, model) -> None:
        variant = model.variants.get("base")
        if variant is None:
            return
        is_active_model = model.model_id == self._state.active_model_id
        active_nid = variant.active_node_id if is_active_model else None

        # Composition
        comp = variant.nodes.get("composition")
        stale = "  *" if (comp and comp.status is NodeStatus.STALE) else ""
        if comp is not None:
            dpg.add_selectable(
                label=f"    Composition{stale}", parent=_MODEL_TREE,
                default_value=(active_nid == "composition"),
                user_data=(model.model_id, "composition"),
                callback=self._on_tree_open_node,
            )

        self._render_model_lab_data_tree(model)

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
            dpg.add_menu_item(label="Delete model...", user_data=model_id,
                              callback=self._on_delete_model_confirm)

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
        self._selected_tree_model_id = mid
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
        self._selected_tree_model_id = None
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
        self._selected_tree_model_id = None
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
    #  Project overview
    # ==================================================================

    def _overview_stats(self, model) -> tuple[int, int, int, int]:
        """Расчёты/attention/running/saved для строки Overview без загрузки модели."""
        if model.loaded and (variant := model.variants.get("base")) is not None:
            nodes = [node for node in variant.nodes.values()
                     if node.kind not in (NodeKind.COMPOSITION, NodeKind.COMPARE)]
            return (
                sum(node.result is not None for node in nodes),
                sum(node.status is NodeStatus.STALE or node.error is not None
                    for node in nodes),
                sum(node.status is NodeStatus.RUNNING for node in nodes),
                0,
            )
        info = proj_svc.calc_summary(model.summary) if model.summary else {}
        return (int(info.get("persisted", 0)), int(info.get("stale", 0)), 0,
                int(info.get("persisted", 0)))

    def _render_project_overview(self, parent) -> None:
        project_id = self._state.active_project_id
        project = self._state.projects.get(project_id) if project_id else None
        models = [self._state.models[mid] for mid in (project.model_ids if project else ())
                  if mid in self._state.models]
        title = project.title if project else "Current project"
        dpg.add_text(f"Project overview: {title}", parent=parent)
        dpg.add_text(
            "Saved results come from the project database; unsaved edits are marked "
            "with * and remain local until Save.",
            parent=parent, wrap=820,
        )
        active_id = self._state.active_model_id
        with dpg.group(horizontal=True, parent=parent):
            dpg.add_button(label="Open active composition", user_data=active_id,
                           callback=self._on_overview_open_composition,
                           enabled=active_id is not None)
            dpg.add_button(label="New flash", user_data=active_id,
                           callback=self._on_overview_new_flash,
                           enabled=active_id is not None)
            dpg.add_button(label="New phase envelope", user_data=active_id,
                           callback=self._on_overview_new_envelope,
                           enabled=active_id is not None)
        dpg.add_separator(parent=parent)
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True,
                       borders_outerV=True, resizable=True, scrollY=True,
                       height=330, freeze_rows=1):
            for label in ("Model", "Field", "Components", "EOS", "Results",
                          "Attention", "Last saved", "Actions"):
                dpg.add_table_column(label=label)
            for model in models:
                results, attention, running, saved = self._overview_stats(model)
                with dpg.table_row():
                    dpg.add_text(model.title + ("  *" if model.dirty else ""))
                    dpg.add_text(model.field_name or "-")
                    dpg.add_text(str(model.n_components))
                    dpg.add_text(model.eos or "-")
                    result_text = str(results)
                    if running:
                        result_text += f" ({running} running)"
                    if saved and not model.loaded:
                        result_text += " saved"
                    dpg.add_text(result_text)
                    attention_text = str(attention) if attention else "-"
                    dpg.add_text(attention_text)
                    saved_at = (model.summary.workspace_saved_at
                                if model.summary is not None else None)
                    dpg.add_text(self._format_saved_at(saved_at) or "-")
                    with dpg.group(horizontal=True):
                        dpg.add_button(label="Open", small=True,
                                       user_data=model.model_id,
                                       callback=self._on_overview_open_model)
                        dpg.add_button(label="Last result", small=True,
                                       user_data=model.model_id,
                                       callback=self._on_overview_open_last_result,
                                       enabled=bool(results))

    def _activate_overview_model(self, model_id: str | None) -> bool:
        if not model_id or model_id not in self._state.models:
            return False
        if model_id != self._state.active_model_id:
            self._state.set_active_model(model_id, notify=False)
            self._restore_workspace(model_id)
        return True

    def _on_overview_open_model(self, sender, app_data, user_data) -> None:
        model_id = str(user_data)
        if not self._activate_overview_model(model_id):
            return
        self._expanded_models.add(model_id)
        self._state.notify(StateChange(StateChangeKind.NAVIGATION))
        self._set_status(f"Active model: '{self._state.models[model_id].title}'.")

    def _on_overview_open_composition(self, sender, app_data, user_data) -> None:
        model_id = str(user_data) if user_data else None
        if self._activate_overview_model(model_id):
            self._state.open_node("composition")

    def _on_overview_new_flash(self, sender, app_data, user_data) -> None:
        model_id = str(user_data) if user_data else None
        if self._activate_overview_model(model_id):
            self._state.new_flash_run()
            self._set_status("New flash tab opened - set P/T and Run.")

    def _on_overview_new_envelope(self, sender, app_data, user_data) -> None:
        model_id = str(user_data) if user_data else None
        if self._activate_overview_model(model_id):
            self._on_new_envelope(None, None, model_id)

    def _on_overview_open_last_result(self, sender, app_data, user_data) -> None:
        model_id = str(user_data)
        if not self._activate_overview_model(model_id):
            return
        variant = self._state.active_variant
        if variant is None:
            return
        results = [node for node in variant.nodes.values()
                   if node.kind not in (NodeKind.COMPOSITION, NodeKind.COMPARE)
                   and node.result is not None]
        if results:
            self._state.open_node(results[-1].node_id)
        else:
            self._set_status("This model has no loaded calculation result.")

    def _overview_is_selected(self) -> bool:
        return bool(self._overview_tab_id is not None and self._tabbar_id is not None
                    and dpg.does_item_exist(self._overview_tab_id)
                    and dpg.get_value(self._tabbar_id) == self._overview_tab_id)

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

        self._tab_ids = {}
        self._tab_content_ids = {}
        self._overview_tab_id = None
        with dpg.tab_bar(parent=_WORKSPACE, reorderable=True,
                         callback=self._on_tab_changed) as tabbar:
            self._tabbar_id = tabbar
            with dpg.tab(label="Project overview") as overview_tab:
                self._overview_tab_id = overview_tab
                self._render_project_overview(overview_tab)
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
        elif self._overview_tab_id is not None:
            dpg.set_value(self._tabbar_id, self._overview_tab_id)

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
        if app_data == self._overview_tab_id:
            self._render_tree()
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
