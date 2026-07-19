"""Headless-smoke реального DearPyGui View без показа окна."""

import shutil
from pathlib import Path

import dearpygui.dearpygui as dpg

from gui.app_state import AppState
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

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"


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


def test_compare_selection_is_scoped_to_model(tmp_path):
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
        labels = [
            dpg.get_item_label(child)
            for children in dpg.get_item_children(_MODEL_TREE).values()
            for child in children
            if dpg.get_item_type(child) == "mvAppItemType::mvSelectable"
        ]
        starred = [label for label in labels if "[*]" in label]
        assert len(starred) == 1
        assert "100 bar / 50 C" in starred[0]
        assert copy_flash == source_flash  # одинаковые локальные id — ключевой случай
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
        # Именно enabled, а не просто наличие: кнопка какое-то время стояла
        # enabled=False, и проверка по подписи этого не замечала.
        paste_btn = _item_by_label(_WORKSPACE, "Paste from Excel")
        assert paste_btn is not None
        assert dpg.get_item_configuration(paste_btn)["enabled"] is True
        assert _has_label(_WORKSPACE, "Copy table")
        assert _has_label(_WORKSPACE, "Bo vs pressure")
        assert _has_label(_WORKSPACE, "vapor_density vs pressure")
        assert len([item for item in dpg.get_item_children(
            app._exp_chart_holder[nid], 1)
            if dpg.get_item_type(item) == "mvAppItemType::mvGroup"]) == 2
        chart_groups = [item for item in dpg.get_item_children(
            app._exp_chart_holder[nid], 1)
            if dpg.get_item_type(item) == "mvAppItemType::mvGroup"]
        app._rebuild_exp_chart_grid(nid)
        assert len([item for item in dpg.get_item_children(
            app._exp_chart_holder[nid], 1)
            if dpg.get_item_type(item) == "mvAppItemType::mvGroup"]) == len(chart_groups)
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
