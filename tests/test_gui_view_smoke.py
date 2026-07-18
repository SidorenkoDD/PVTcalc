"""Headless-smoke реального DearPyGui View без показа окна."""

from pathlib import Path

import dearpygui.dearpygui as dpg

from gui.app_state import AppState
from gui.services.model_repository import ModelRepository
from gui.session import SessionState
from gui.view.app import _PRIMARY, _PROJECTS_CONTENT, PVTcalcApp

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
        assert dpg.does_item_exist(_PRIMARY)
        assert dpg.does_item_exist(_PROJECTS_CONTENT)
        assert state.models
    finally:
        dpg.destroy_context()
