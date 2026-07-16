"""
Тесты undo/redo (снимки состояния варианта). Без DearPyGui — идут в CI.
"""

from pathlib import Path

import pytest

from gui.app_state import AppState, NodeStatus
from gui.services.model_repository import ModelRepository

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"


@pytest.fixture
def state():
    st = AppState(ModelRepository(db_path=str(MODELS_JSON)))
    st.refresh_model_list()
    st.open_model("KRSNL_PVTSIM")
    st.open_node("composition")
    return st


def test_undo_delete_restores_node(state):
    a = state.new_flash_run(50, 20)
    assert a in state.active_variant.nodes
    state.delete_node(a)
    assert a not in state.active_variant.nodes
    state.undo()
    assert a in state.active_variant.nodes  # узел вернулся


def test_redo_after_undo(state):
    a = state.new_flash_run(50, 20)
    state.delete_node(a)
    state.undo()
    assert a in state.active_variant.nodes
    state.redo()
    assert a not in state.active_variant.nodes  # удаление повторилось


def test_undo_composition_edit(state):
    before = state.active_composition.composition["C1"]
    state.edit_zi("C1", 0.999)
    assert state.active_composition.composition["C1"] == 0.999
    state.undo()
    assert state.active_composition.composition["C1"] == before


def test_undo_eos_change(state):
    assert state.active_composition.eos_name.value == "PREOS"
    state.set_composition_eos("SRKEOS")
    assert state.active_composition.eos_name.value == "SRKEOS"
    state.undo()
    assert state.active_composition.eos_name.value == "PREOS"


def test_undo_rename(state):
    a = state.new_flash_run(50, 20)
    state.rename_node(a, "bubble")
    assert state.node_by_id(a).params["label"] == "bubble"
    state.undo()
    assert "label" not in state.node_by_id(a).params


def test_new_edit_clears_redo(state):
    a = state.new_flash_run(50, 20)
    state.delete_node(a)
    state.undo()
    assert state.can_redo()
    state.new_flash_run(100, 30)  # новое действие очищает redo
    assert not state.can_redo()


def test_can_undo_redo_flags(state):
    assert not state.can_undo()
    state.new_flash_run(50, 20)
    assert state.can_undo()
    assert not state.can_redo()
    state.undo()
    assert state.can_redo()


def test_undo_is_per_variant(state):
    # у только что открытой модели история пуста
    state.new_flash_run(50, 20)
    assert state.can_undo()
    # другая модель — своя (пустая) история
    state.open_model("PRRZLM_MDT_TEST")
    assert not state.can_undo()
