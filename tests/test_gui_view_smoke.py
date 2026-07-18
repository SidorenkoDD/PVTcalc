"""Headless-smoke реального DearPyGui View без показа окна."""

from pathlib import Path

import dearpygui.dearpygui as dpg

from gui.app_state import AppState
from gui.services.model_repository import ModelRepository
from gui.session import SessionState
from gui.view.app import _PRIMARY, _PROJECTS_CONTENT, _WORKSPACE, PVTcalcApp
from gui.view.contracts import ViewHost

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"


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

        # Реально собрать вкладки каждого вынесенного view-модуля. Расчёты не
        # запускаются: smoke проверяет DPG wiring/callback construction.
        model_id = next(iter(state.models))
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
