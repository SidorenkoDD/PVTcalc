"""Headless-smoke реального DearPyGui View без показа окна."""

from pathlib import Path

import dearpygui.dearpygui as dpg

from gui.app_state import AppState
from gui.services.model_repository import ModelRepository
from gui.session import SessionState
from gui.view.app import _PRIMARY, _PROJECTS_CONTENT, _WORKSPACE, PVTcalcApp
from gui.view.contracts import ViewHost

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
    pending: list[int | str] = [root]
    while pending:
        item = pending.pop()
        for children in dpg.get_item_children(item).values():
            pending.extend(children)
            if any(dpg.get_item_label(child) == label for child in children):
                return True
    return False


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
        app._on_duplicate_model_confirm(None, None, model_id)
        assert app._duplicate_win is not None
        assert _has_label(app._duplicate_win, "Model name")
        app._close_duplicate_modal()

        # Реально собрать вкладки каждого вынесенного view-модуля. Расчёты не
        # запускаются: smoke проверяет DPG wiring/callback construction.
        app._open_project(model_id)
        assert dpg.does_item_exist(_WORKSPACE)
        state.open_node("composition")
        assert "composition" in app._tab_ids

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
        assert _has_label(_WORKSPACE, "Bo vs pressure")
        assert _has_label(_WORKSPACE, "vapor_density vs pressure")
        assert len([item for item in dpg.get_item_children(
            app._exp_chart_holder[nid], 1)
            if dpg.get_item_type(item) == "mvAppItemType::mvGroup"]) == 2
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
