"""Headless-smoke реального DearPyGui View без показа окна."""

import json
import shutil
from pathlib import Path

import dearpygui.dearpygui as dpg

from gui.app_state import AppState
from gui.services import lab_data_service as lab_svc
from gui.services import project_service as proj_svc
from gui.services.model_repository import ModelRepository
from gui.session import SessionState
from gui.view.app import (
    _PRIMARY,
    _PROJECTS_CONTENT,
    _STATUS_BAR,
    _WORKSPACE,
    PVTcalcApp,
)
from gui.view.contracts import ViewHost
from gui.view.workspace_view import _MODEL_TREE

MODELS_JSON = Path(__file__).resolve().parent / "fixtures" / "models.json"


def _text_values(root: int | str) -> list[str]:
    """Собирает значения DPG text-виджетов под контейнером ``root``."""
    values: list[str] = []
    pending: list[int | str] = [root]
    while pending:
        item = pending.pop()
        for children in dpg.get_item_children(item).values():
            pending.extend(children)
            for child in children:
                if dpg.get_item_type(child) == "mvAppItemType::mvText":
                    values.append(str(dpg.get_value(child)))
    return values


def _has_label(root: int | str, label: str) -> bool:
    return _item_by_label(root, label) is not None


def _item_by_label(root: int | str, label: str) -> int | None:
    """Первый виджет под ``root`` с подписью ``label`` (или None).

    В отличие от `_has_label` возвращает сам id — чтобы проверять не только
    наличие виджета, но и его конфигурацию (например, enabled у кнопки).
    """
    pending: list[int | str] = [root]
    while pending:
        item = pending.pop()
        if not dpg.does_item_exist(item):
            continue
        for children in dpg.get_item_children(item).values():
            pending.extend(children)
            for child in children:
                if dpg.does_item_exist(child) and dpg.get_item_label(child) == label:
                    return child
    return None


def test_projects_renders_corrupt_store_diagnostic(tmp_path):
    path = tmp_path / "models.json"
    path.write_text('{"broken":', encoding="utf-8")
    state = AppState(ModelRepository(str(path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        texts = _text_values(_PROJECTS_CONTENT)
        assert any("Models database could not be read" in text for text in texts)
        assert any("строка 1, столбец 11" in text for text in texts)
        assert _has_label(_PROJECTS_CONTENT, "Refresh models")
    finally:
        dpg.destroy_context()


def test_dpg_context_builds_and_renders_projects():
    state = AppState(ModelRepository(str(MODELS_JSON)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=800, height=600)
        app._build_menu()
        app._build_layout()
        app._build_shortcuts()
        state.refresh_model_list()
        assert isinstance(app, ViewHost)
        assert app._view_context.state is state
        assert app._view_context.jobs is app._jobs
        assert dpg.does_item_exist(_PRIMARY)
        assert dpg.does_item_exist(_PROJECTS_CONTENT)
        assert state.models
        model_id = next(iter(state.models))

        # Реально собрать вкладки каждого вынесенного view-модуля. Расчёты не
        # запускаются: smoke проверяет DPG wiring/callback construction.
        app._open_project(model_id)
        assert dpg.does_item_exist(_WORKSPACE)
        assert any(text.startswith("Models in project")
                   for text in _text_values(_MODEL_TREE))
        assert any(text.startswith("Shortcuts: Ctrl+S save")
                   for text in _text_values(_MODEL_TREE))
        root_labels = [
            dpg.get_item_label(child)
            for children in dpg.get_item_children(_MODEL_TREE).values()
            for child in children
            if dpg.get_item_type(child) == "mvAppItemType::mvSelectable"
        ]
        assert any(state.models[model_id].title in label for label in root_labels)
        assert any(label.startswith("[A] ") and state.models[model_id].title in label
                   for label in root_labels)
        state.open_node("composition")
        assert "composition" in app._tab_ids
        app._on_duplicate_model_confirm(None, None, model_id)
        assert app._duplicate_win is not None
        assert _has_label(app._duplicate_win, "Model name")
        app._close_duplicate_modal()

        flash_id = state.new_flash_run(100.0, 50.0)
        experiment_id = state.new_experiment("dle", {"pressures": [200.0, 100.0]})
        envelope_id = state.new_envelope({"method": "ssm"})
        assert {flash_id, experiment_id, envelope_id}.issubset(app._tab_ids)

        app._on_open_settings()
        assert dpg.does_item_exist(app._settings_win)
        assert app._settings_ids
        assert all(dpg.get_item_configuration(wid)["readonly"]
                   for wid in app._settings_ids.values())
    finally:
        dpg.destroy_context()


def test_workspace_tree_groups_duplicate_models_by_project(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        labels = [
            dpg.get_item_label(child)
            for children in dpg.get_item_children(_MODEL_TREE).values()
            for child in children
            if dpg.get_item_type(child) == "mvAppItemType::mvSelectable"
        ]
        assert any(state.models["KRSNL_PVTSIM"].title in label for label in labels)
        assert any(state.models[copy_id].title in label for label in labels)
    finally:
        dpg.destroy_context()


def test_project_overview_lists_models_and_opens_selected_model(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=1000, height=700)
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")

        assert _has_label(_WORKSPACE, "Project overview")
        texts = _text_values(_WORKSPACE)
        assert any("Project overview:" in text for text in texts)
        assert state.models[copy_id].title in texts
        assert _has_label(_WORKSPACE, "Open active composition")
        assert _has_label(_WORKSPACE, "New flash")

        app._on_overview_open_model(None, None, copy_id)

        assert state.active_model_id == copy_id
        assert copy_id in app._expanded_models
    finally:
        dpg.destroy_context()


def test_tree_search_groups_dle_runs_by_source_model(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        first_id = state.new_experiment("dle", {"T_c": 50.0, "P_res": 100.0})
        assert first_id is not None
        app._on_model_row(None, None, copy_id)
        second_id = state.new_experiment("dle", {"T_c": 60.0, "P_res": 120.0})
        assert second_id is not None

        app._on_tree_query(None, "DLE")

        labels: list[str] = []
        pending: list[int | str] = [_MODEL_TREE]
        while pending:
            for children in dpg.get_item_children(pending.pop()).values():
                pending.extend(children)
                labels.extend(
                    dpg.get_item_label(child) for child in children
                    if dpg.get_item_type(child) == "mvAppItemType::mvSelectable"
                )
        assert any(state.models["KRSNL_PVTSIM"].title in label for label in labels)
        assert any(state.models[copy_id].title in label for label in labels)
        assert sum("DLE —" in label for label in labels) == 2

        app._on_tree_open_node(None, None, (copy_id, second_id))

        assert state.active_model_id == copy_id
        assert state.active_variant is not None
        assert state.active_variant.active_node_id == second_id
    finally:
        dpg.destroy_context()


def test_search_selects_models_and_opens_cross_model_compare(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        params = {"pressures": [300, 200], "T_c": 50.0, "P_res": 300.0}
        first_id = state.new_experiment("dle", params)
        assert first_id is not None
        app._on_model_row(None, None, copy_id)
        second_id = state.new_experiment("dle", params)
        assert second_id is not None

        app._on_tree_query(None, "DLE")
        app._on_search_model_selected(None, True, "KRSNL_PVTSIM")
        app._on_search_model_selected(None, True, copy_id)
        assert app._search_selected_model_ids == {"KRSNL_PVTSIM", copy_id}

        app._on_compare_search_models(None, None)

        compare = state.node_by_id("compare")
        assert compare is not None
        assert [member["model_id"] for member in compare.params["member_refs"]] == [
            "KRSNL_PVTSIM", copy_id,
        ]
        assert len(app._compare_selection) == 2
    finally:
        dpg.destroy_context()


def test_project_lab_tree_autosaves_only_nonempty_manual_dataset(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")

        app._on_new_lab_dataset(None, None, ("dle", "project", None))
        assert app._lab_catalog_editor is not None
        assert app._lab_catalog_editor["title"] == "DLE LAB-1"
        assert _item_by_label(app._lab_catalog_modal, "Save") is None
        assert _has_label(app._lab_catalog_modal, "Close")
        assert _has_label(app._lab_catalog_modal, "T res, C (optional)")
        assert lab_svc.list_datasets(str(db_path), "KRSNL_PVTSIM",
                                     experiment_kind="dle") == []
        app._on_catalog_lab_add_row()
        columns = app._lab_catalog_editor["columns"]
        for index in range(len(columns)):
            app._on_catalog_lab_cell(None, str(300 - index), (0, index))

        datasets = lab_svc.list_datasets(str(db_path), "KRSNL_PVTSIM",
                                         experiment_kind="dle")
        assert len(datasets) == 1
        assert datasets[0]["rows"][0] == [float(300 - i)
                                           for i in range(len(columns))]
        assert app._lab_catalog_editor is not None

        # Del в открытом редакторе вызывает то же подтверждение, что и меню
        # дерева: отдельной кнопки удаления в редакторе больше нет.
        app._on_key_delete(None, None)
        assert any("Linked experiments will lose this source" in text
                   for text in _text_values(app._modals[-1]))
    finally:
        dpg.destroy_context()


def test_lab_dataset_single_click_selects_and_double_click_opens_editor(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    dataset = lab_svc.create_dataset(
        str(db_path), "KRSNL_PVTSIM", title="DLE LAB-1", experiment_kind="dle",
        columns=["pressure"], rows=[[300.0]],
    )
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        ref = (dataset["dataset_id"], "project", None)

        app._on_project_lab_dataset(None, True, ref)
        assert app._selected_tree_lab_dataset == ref
        assert app._lab_catalog_modal is None

        app._on_project_lab_dataset_edit(None, None, ref)
        assert app._lab_catalog_modal is not None
        assert dpg.does_item_exist(app._lab_catalog_modal)

        # Closing via the regular action must clear both the window id and
        # editor state, otherwise a later double-click can be swallowed.
        app._on_catalog_lab_cancel()
        assert app._lab_catalog_modal is None
        assert app._lab_catalog_editor is None
        app._on_project_lab_dataset_edit(None, None, ref)
        assert app._lab_catalog_modal is not None
        assert dpg.does_item_exist(app._lab_catalog_modal)
    finally:
        dpg.destroy_context()


def test_lab_data_tree_groups_datasets_by_experiment_kind(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    common = {"columns": ["pressure"], "rows": [[300.0]]}
    lab_svc.create_dataset(str(db_path), "KRSNL_PVTSIM", title="DLE A",
                           experiment_kind="dle", **common)
    lab_svc.create_dataset(str(db_path), "KRSNL_PVTSIM", title="DLE B",
                           experiment_kind="dle", **common)
    lab_svc.create_dataset(str(db_path), "KRSNL_PVTSIM", title="CCE A",
                           experiment_kind="cce", **common)
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        app._expanded_cats.update({"KRSNL_PVTSIM:lab", "KRSNL_PVTSIM:lab:dle"})
        app._render_tree()

        labels = [
            dpg.get_item_label(child)
            for children in dpg.get_item_children(_MODEL_TREE).values()
            for child in children
            if dpg.get_item_type(child) == "mvAppItemType::mvSelectable"
        ]
        assert "    v DLE (2)" in labels
        assert "    > CCE (1)" in labels
        assert any("DLE A" in label for label in labels)
        assert not any("?" in label for label in labels if "DLE A" in label)
    finally:
        dpg.destroy_context()


def test_workspace_tree_keeps_multiple_models_expanded(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")

        # Первая модель раскрыта при входе; один клик по второй одновременно
        # активирует и раскрывает её, не схлопывая первую.
        app._on_model_row(None, None, copy_id)
        assert {"KRSNL_PVTSIM", copy_id} <= app._expanded_models
        root_labels = [
            dpg.get_item_label(child)
            for children in dpg.get_item_children(_MODEL_TREE).values()
            for child in children
            if dpg.get_item_type(child) == "mvAppItemType::mvSelectable"
        ]
        assert sum("Composition" in label for label in root_labels) == 2

        # Повторный одинарный клик сразу меняет стрелку и схлопывает активную
        # модель, при этом исходная остаётся раскрытой.
        app._on_model_row(None, None, copy_id)
        assert copy_id not in app._expanded_models
        assert "KRSNL_PVTSIM" in app._expanded_models
    finally:
        dpg.destroy_context()


def test_tree_category_click_activates_its_model(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        app._on_model_row(None, None, copy_id)  # загрузить и раскрыть копию

        # Оставляем обе модели раскрытыми, но активируем исходную перед кликом
        # по категории копии — это повторяет реальное IDE-переключение.
        app._state.set_active_model("KRSNL_PVTSIM", notify=False)
        app._restore_workspace("KRSNL_PVTSIM")
        app._render_tree()
        app._on_cat_toggle(None, None, f"{copy_id}:flash")

        assert state.active_model_id == copy_id
        assert f"{copy_id}:flash" in app._expanded_cats
        assert dpg.get_value(_STATUS_BAR).startswith("Active model:")
    finally:
        dpg.destroy_context()


def test_model_tree_context_save_targets_requested_model(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        model = state.models["KRSNL_PVTSIM"]
        model.dirty = True
        app._render_tree()

        app._on_save_tree_model(None, None, model.model_id)

        assert model.dirty is False
        assert dpg.get_value(_STATUS_BAR) == "Model 'KRSNL_PVTSIM' saved."
    finally:
        dpg.destroy_context()


def test_project_delete_dialog_warns_about_dirty_loaded_models(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=800, height=600)
        app._build_layout()
        state.refresh_model_list()
        state.enter_model("KRSNL_PVTSIM")
        state.models["KRSNL_PVTSIM"].dirty = True
        project_id = state.models["KRSNL_PVTSIM"].project_id

        app._on_delete_project_confirm(None, None, project_id)

        modal = app._modals[-1]
        assert any("1 model(s) have unsaved changes" in text
                   for text in _text_values(modal))
        dpg.delete_item(modal)
    finally:
        dpg.destroy_context()


def test_back_to_projects_confirms_and_saves_full_workspace(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=800, height=600)
        app._build_layout()
        state.refresh_model_list()
        state.enter_model("KRSNL_PVTSIM")
        node_id = state.new_envelope({"method": "ssm", "t_min_c": 10.0})
        assert node_id is not None
        state.set_node_result(node_id, {"points": [[10.0, 42.0]]})

        app._on_back_to_projects()

        assert state.current_screen == "workspace"
        modal = app._modals[-1]
        assert any("Save 1 changed model(s)" in text for text in _text_values(modal))
        app._save_then_continue(modal, state.show_projects)

        assert state.current_screen == "projects"
        raw = json.loads(db_path.read_text(encoding="utf-8"))
        assert raw["KRSNL_PVTSIM"]["workspace"]["snapshot"]["nodes"]
        assert raw["KRSNL_PVTSIM"]["results"][0]["has_result"] is True
    finally:
        dpg.destroy_context()


def test_discard_restores_previous_saved_workspace_not_empty_model(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    # Сначала создаём именно постоянный расчёт, который пользователь откроет
    # в следующей сессии.
    seed = AppState(ModelRepository(str(db_path)))
    seed.refresh_model_list()
    seed.enter_model("KRSNL_PVTSIM")
    saved_id = seed.new_envelope({"method": "ssm", "t_min_c": 10.0})
    assert saved_id == "env_1"
    seed.set_node_result(saved_id, {"points": [[10.0, 42.0]]})
    seed.save_model()

    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=800, height=600)
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        assert state.node_by_id(saved_id).result == {"points": [[10.0, 42.0]]}

        unsaved_id = state.new_envelope({"method": "grid", "grid_t_points": 5})
        assert unsaved_id == "env_2"
        state.set_node_result(unsaved_id, {"points": [[20.0, 99.0]]})
        app._on_back_to_projects()
        modal = app._modals[-1]
        app._discard_then_continue(modal, state.show_projects)

        app._open_project("KRSNL_PVTSIM")
        restored = state.active_variant.nodes
        assert restored[saved_id].result == {"points": [[10.0, 42.0]]}
        assert unsaved_id not in restored
    finally:
        dpg.destroy_context()


def test_delete_model_confirmation_warns_and_removes_only_that_model(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState(workspaces={
        "KRSNL_PVTSIM": {"legacy": "will be removed"},
    }))
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=800, height=600)
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        saved_id = state.new_envelope({"method": "ssm", "t_min_c": 10.0})
        assert saved_id is not None
        state.set_node_result(saved_id, {"points": [[10.0, 42.0]]})
        state.save_model()
        state.new_flash_run()  # несохранённое изменение тоже должно быть уничтожено

        app._on_delete_model_confirm(None, None, "KRSNL_PVTSIM")

        modal = app._modals[-1]
        texts = _text_values(modal)
        assert any("all saved calculation results" in text for text in texts)
        assert any("unsaved changes will be discarded" in text for text in texts)
        app._on_delete_model_do(None, None, ("KRSNL_PVTSIM", modal))

        raw = json.loads(db_path.read_text(encoding="utf-8"))
        assert "KRSNL_PVTSIM" not in raw
        assert copy_id in raw
        assert "KRSNL_PVTSIM" not in state.models
        assert state.active_model_id == copy_id
        assert state.current_screen == "workspace"
        assert app._session.workspaces == {}
    finally:
        dpg.destroy_context()


def test_delete_key_for_selected_model_opens_confirmation(tmp_path, monkeypatch):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=800, height=600)
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        app._on_model_row(None, None, "KRSNL_PVTSIM")
        root_id = next(
            child for children in dpg.get_item_children(_MODEL_TREE).values()
            for child in children
            if (dpg.get_item_type(child) == "mvAppItemType::mvSelectable"
                and state.models["KRSNL_PVTSIM"].title in dpg.get_item_label(child))
        )
        monkeypatch.setattr(dpg, "is_item_hovered", lambda _item: False)
        monkeypatch.setattr(dpg, "get_focused_item", lambda: root_id)

        app._on_key_delete(None, None)

        assert "KRSNL_PVTSIM" in state.models  # Del сначала требует подтверждения
        assert any("Delete model" in text for text in _text_values(app._modals[-1]))
    finally:
        dpg.destroy_context()


def test_global_edit_shortcuts_do_not_run_behind_modal(monkeypatch):
    state = AppState(ModelRepository(str(MODELS_JSON)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        with dpg.window(label="test modal", modal=True) as modal:
            app._track_modal(modal)
            calls: list[str] = []
            monkeypatch.setattr(app, "_on_save_model", lambda: calls.append("save"))
            monkeypatch.setattr(app, "_do_undo", lambda: calls.append("undo"))
            monkeypatch.setattr(app, "_do_redo", lambda: calls.append("redo"))
            monkeypatch.setattr(dpg, "is_key_down", lambda _key: True)

            app._on_key_s(None, None)
            app._on_key_z(None, None)
            app._on_key_y(None, None)

            assert calls == []
    finally:
        dpg.destroy_context()


def test_system_exit_request_uses_save_confirmation(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=800, height=600)
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        assert state.new_flash_run() is not None

        app._on_viewport_close_request()

        modal = app._modals[-1]
        assert any("before closing the application" in text
                   for text in _text_values(modal))
        assert _has_label(modal, "Save all")
        assert _has_label(modal, "Discard changes")
        assert _has_label(modal, "Cancel")
    finally:
        dpg.destroy_context()


def test_exit_request_dismisses_lab_editor_then_shows_save_confirmation(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=800, height=600)
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        assert state.new_flash_run() is not None
        app._on_new_lab_dataset(None, None, ("dle", "project", None))
        assert app._lab_catalog_modal is not None

        app._on_viewport_close_request()

        assert app._lab_catalog_modal is None
        modal = app._modals[-1]
        assert any("before closing the application" in text
                   for text in _text_values(modal))
    finally:
        dpg.destroy_context()


def test_compare_selection_distinguishes_same_local_id_across_models(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        app._expanded_cats.update({"KRSNL_PVTSIM:flash", f"{copy_id}:flash"})

        source_flash = state.new_flash_run(100.0, 50.0)
        assert source_flash == "flash_1"
        app._on_model_row(None, None, copy_id)
        copy_flash = state.new_flash_run(100.0, 50.0)
        assert copy_flash == "flash_1"

        app._toggle_compare("KRSNL_PVTSIM", source_flash)
        app._toggle_compare(copy_id, copy_flash)
        labels = [
            dpg.get_item_label(child)
            for children in dpg.get_item_children(_MODEL_TREE).values()
            for child in children
            if dpg.get_item_type(child) == "mvAppItemType::mvSelectable"
        ]
        starred = [label for label in labels if "[*]" in label]
        assert len(starred) == 2
        assert all("100 bar / 50 C" in label for label in starred)
        assert len(app._compare_selection) == 2
        assert {ref.model_id for ref in app._compare_selection} == {
            "KRSNL_PVTSIM", copy_id,
        }
        assert copy_flash == source_flash  # одинаковые локальные id — ключевой случай
    finally:
        dpg.destroy_context()


def test_same_experiment_id_opens_after_switching_between_duplicate_models(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")

        first_id = state.new_experiment(
            "dle", {"pressures": [400, 200], "T_c": 100.0, "P_res": 400.0},
        )
        app._on_model_row(None, None, copy_id)
        second_id = state.new_experiment(
            "dle", {"pressures": [400, 200], "T_c": 100.0, "P_res": 400.0},
        )
        assert first_id == second_id == "exp_1"

        app._on_tree_open_node(None, None, ("KRSNL_PVTSIM", first_id))
        assert state.active_model_id == "KRSNL_PVTSIM"
        assert state.active_variant.active_node_id == first_id
        assert first_id in app._tab_ids
        assert dpg.does_item_exist(app._tab_ids[first_id])

        app._on_tree_open_node(None, None, (copy_id, second_id))
        assert state.active_model_id == copy_id
        assert state.active_variant.active_node_id == second_id
        assert second_id in app._tab_ids
        assert dpg.does_item_exist(app._tab_ids[second_id])
    finally:
        dpg.destroy_context()


def test_compare_renders_same_experiment_from_two_models(tmp_path):
    db_path = tmp_path / "models.json"
    shutil.copyfile(MODELS_JSON, db_path)
    copy_id = proj_svc.duplicate_model(str(db_path), "KRSNL_PVTSIM")
    state = AppState(ModelRepository(str(db_path)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        params = {"pressures": [300, 200], "T_c": 100.0, "P_res": 400.0}
        result = {
            "kind": "dle", "x": "pressure",
            "columns": ["pressure", "Bo", "Rs"],
            "rows": [[300.0, 1.3, 120.0], [200.0, 1.2, 100.0]],
            "plot_all": ["Bo", "Rs"], "charts": ["Bo", "Rs"],
            "main_columns": ["pressure", "Bo", "Rs"], "stages": [],
        }
        first_id = state.new_experiment("dle", params)
        first_ref = state.node_ref(first_id)
        state.set_node_result(first_ref, result)

        app._on_model_row(None, None, copy_id)
        second_id = state.new_experiment("dle", params)
        second_ref = state.node_ref(second_id)
        tuned = {**result, "rows": [[300.0, 1.35, 123.0],
                                     [200.0, 1.18, 98.0]]}
        state.set_node_result(second_ref, tuned)

        app._toggle_compare("KRSNL_PVTSIM", first_id)
        app._toggle_compare(copy_id, second_id)
        app._on_open_compare(None, None, copy_id)

        compare = state.node_by_id("compare")
        assert compare is not None
        assert [item["model_id"] for item in compare.params["member_refs"]] == [
            "KRSNL_PVTSIM", copy_id,
        ]
        resolved = app._resolve_compare_members(compare)
        assert len(resolved) == 2
        labels = [app._exp_compare_label(ref, node) for ref, node in resolved]
        assert state.models["KRSNL_PVTSIM"].title in labels[0]
        assert state.models[copy_id].title in labels[1]
        assert "compare" in app._tab_ids
        assert dpg.does_item_exist(app._tab_ids["compare"])
    finally:
        dpg.destroy_context()


def test_node_status_update_keeps_other_workspace_tabs_intact():
    state = AppState(ModelRepository(str(MODELS_JSON)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        state.enter_model(next(iter(state.models)))
        state.open_node("composition")
        flash_id = state.new_flash_run(100.0, 20.0)
        assert flash_id is not None

        tabbar = app._tabbar_id
        composition_page = app._tab_content_ids["composition"]
        flash_page = app._tab_content_ids[flash_id]
        state.set_node_running(flash_id)

        assert app._tabbar_id == tabbar
        assert app._tab_content_ids["composition"] == composition_page
        assert app._tab_content_ids[flash_id] == flash_page
        assert state.node_by_id(flash_id).status.name == "RUNNING"
    finally:
        dpg.destroy_context()


def test_experiment_lab_data_and_chart_grid_render():
    state = AppState(ModelRepository(str(MODELS_JSON)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        state.enter_model(next(iter(state.models)))
        state.open_node("composition")
        nid = state.new_experiment("dle", {"pressures": [400.0, 200.0]})
        node = state.node_by_id(nid)
        node.result = {
            "columns": ["pressure", "Bo", "Rs", "liquid_density",
                        "vapor_density"],
            "rows": [[400.0, 1.3, 90.0, 700.0, 20.0],
                     [200.0, 1.1, 70.0, 710.0, 18.0]],
            "x": "pressure",
            "charts": ["Bo", "Rs", "liquid_density", "vapor_density"],
            "plot_all": ["Bo", "Rs", "liquid_density", "vapor_density"],
            "main_columns": ["pressure", "Bo", "Rs"],
            "stages": [],
        }
        state.add_lab_data_row(nid, ["pressure", "Bo", "Rs", "liquid_density",
                                     "vapor_density"])
        state.set_lab_data_value(nid, 0, 0, 350.0)
        state.set_lab_data_value(nid, 0, 1, 1.2)
        app._render_node_content(nid)

        assert _has_label(_WORKSPACE, "Lab Data (measured)")
        assert _item_by_label(_WORKSPACE, "Paste from Excel") is None
        assert any("Legacy local Lab Data is read only" in text
                   for text in _text_values(_WORKSPACE))
        assert _has_label(_WORKSPACE, "Copy table")
        assert _has_label(_WORKSPACE, "Bo vs pressure")
        assert _has_label(_WORKSPACE, "vapor_density vs pressure")
        chart_groups = [item for item in dpg.get_item_children(
            app._exp_chart_holder[nid], 1)
            if dpg.get_item_type(item) == "mvAppItemType::mvGroup"]
        assert len(chart_groups) == 4 // app._chart_grid_columns()
        app._rebuild_exp_chart_grid(nid)
        assert len([item for item in dpg.get_item_children(
            app._exp_chart_holder[nid], 1)
            if dpg.get_item_type(item) == "mvAppItemType::mvGroup"]) == len(chart_groups)
    finally:
        dpg.destroy_context()


def test_chart_layout_expands_on_wide_viewport_and_rebuilds_after_resize(monkeypatch):
    state = AppState(ModelRepository(str(MODELS_JSON)))
    app = PVTcalcApp(state, SessionState(window_width=1920, window_height=1080))
    dpg.create_context()
    try:
        dpg.create_viewport(title="test", width=1920, height=1080)
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        nid = state.new_experiment("dle", {"pressures": [400.0, 200.0]})

        assert app._chart_grid_columns() == 2
        assert app._chart_card_width(2) > 440
        assert app._envelope_plot_height() > 520

        calls: list[str] = []
        callbacks = []
        monkeypatch.setattr(app, "_render_node_content", calls.append)
        monkeypatch.setattr(dpg, "get_frame_count", lambda: 100)
        monkeypatch.setattr(dpg, "set_frame_callback",
                            lambda _frame, callback: callbacks.append(callback))
        app._on_viewport_resize()
        callbacks[-1]()

        assert calls == [nid]
    finally:
        dpg.destroy_context()


def test_lab_paste_empty_clipboard_is_safe_and_column_is_vertical(monkeypatch):
    state = AppState(ModelRepository(str(MODELS_JSON)))
    app = PVTcalcApp(state, SessionState())
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        state.enter_model(next(iter(state.models)))
        state.open_node("composition")
        nid = state.new_experiment("dle", {"pressures": [400.0, 200.0]})
        columns = ["pressure", "Bo", "Rs", "liquid_density", "vapor_density"]
        for _ in range(4):
            state.add_lab_data_row(nid, columns)
        app._render_node_content(nid)

        monkeypatch.setattr(dpg, "get_clipboard_text", lambda: "")
        app._on_lab_paste(None, None, nid)
        assert dpg.get_value(_STATUS_BAR).startswith("Clipboard is empty")

        monkeypatch.setattr(dpg, "get_clipboard_text",
                            lambda: "350\n300\n250\n200")
        app._on_lab_paste(None, None, nid)
        rows = state.node_by_id(nid).params["lab_data"]["rows"]
        assert [row[0] for row in rows] == [350.0, 300.0, 250.0, 200.0]
        assert len(app._lab_cell_ids) == 20
        assert app._lab_navigation_registry_id is not None

        current_id = app._lab_cell_ids[(nid, 1, 0)]
        next_id = app._lab_cell_ids[(nid, 2, 0)]
        monkeypatch.setattr(dpg, "get_focused_item", lambda: current_id)
        focused: list[int] = []
        monkeypatch.setattr(dpg, "focus_item", focused.append)
        app._lab_active_cell = (nid, 1, 0)
        app._on_lab_arrow_key(None, None, (1, 0))
        assert focused == [next_id]
        assert app._lab_active_cell == (nid, 2, 0)
        second_id = app._lab_cell_ids[(nid, 3, 0)]
        monkeypatch.setattr(dpg, "get_focused_item", lambda: next_id)
        app._on_lab_arrow_key(None, None, (1, 0))
        assert focused == [next_id, second_id]
        assert app._lab_active_cell == (nid, 3, 0)

        # Native Ctrl+V into an InputText first arrives as one app_data value;
        # it must be expanded into rows instead of remaining in that cell.
        state.clear_lab_data(nid)
        for _ in range(4):
            state.add_lab_data_row(nid, columns)
        app._rebuild_lab_data_table(nid)
        cell_id = app._lab_cell_ids[(nid, 0, 0)]
        app._on_lab_cell(cell_id, "310,280,250,220",
                         (nid, 0, 0))
        rows = state.node_by_id(nid).params["lab_data"]["rows"]
        assert [row[0] for row in rows] == [310.0, 280.0, 250.0, 220.0]

        # Ctrl+V is also handled before InputText can collapse the clipboard
        # range into one widget value.
        state.clear_lab_data(nid)
        for _ in range(4):
            state.add_lab_data_row(nid, columns)
        app._rebuild_lab_data_table(nid)
        current_id = app._lab_cell_ids[(nid, 0, 0)]
        monkeypatch.setattr(dpg, "get_focused_item", lambda: current_id)
        monkeypatch.setattr(dpg, "get_value", lambda item: "")
        monkeypatch.setattr(dpg, "is_key_down",
                            lambda key: key == dpg.mvKey_ModCtrl)
        monkeypatch.setattr(dpg, "get_clipboard_text",
                            lambda: "315,285,255,225")
        app._on_lab_ctrl_v(None, None, None)
        rows = state.node_by_id(nid).params["lab_data"]["rows"]
        assert [row[0] for row in rows] == [315.0, 285.0, 255.0, 225.0]
    finally:
        dpg.destroy_context()


def test_projects_open_survives_malformed_saved_workspace():
    state = AppState(ModelRepository(str(MODELS_JSON)))
    session = SessionState(workspaces={
        "KRSNL_PVTSIM": {
            "nodes": [None, {"kind": "UNKNOWN"}],
            "sequences": {"flash": "bad"},
            "open_tabs": [None],
        },
    })
    app = PVTcalcApp(state, session)
    dpg.create_context()
    try:
        app._build_layout()
        state.refresh_model_list()
        app._open_project("KRSNL_PVTSIM")
        assert state.current_screen == "workspace"
        assert state.active_model_id == "KRSNL_PVTSIM"
    finally:
        dpg.destroy_context()
