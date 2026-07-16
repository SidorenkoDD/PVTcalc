"""
Тесты Фазы 3: сервис флэша (`gui.services.flash_service`) и команды/статусы
узла «Флэш» в `AppState`. Без DearPyGui и без потоков — идут в CI.
"""

from pathlib import Path

import pytest

from calc_core.VLE.FlashResult import FlashResult
from gui.app_state import AppState, NodeKind, NodeStatus
from gui.services import flash_service
from gui.services.model_repository import ModelRepository

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"


@pytest.fixture
def repo() -> ModelRepository:
    return ModelRepository(db_path=str(MODELS_JSON))


@pytest.fixture
def state(repo):
    st = AppState(repo)
    st.refresh_model_list()
    st.open_model("KRSNL_PVTSIM")
    return st


# --- сервис ----------------------------------------------------------------

def test_run_flash_two_phase(repo):
    result = flash_service.run_flash(repo.load_composition("KRSNL_PVTSIM"), 50, 20)
    assert isinstance(result, FlashResult)
    assert result.is_two_phase is True
    assert 0.0 < result.vapor.mole_fraction < 1.0
    assert result.vapor.mole_fraction + result.liquid.mole_fraction == pytest.approx(1.0)
    # свойства обеих фаз посчитаны
    assert result.vapor.properties.get("density") is not None
    assert result.liquid.properties.get("density") is not None


def test_run_flash_single_phase(repo):
    result = flash_service.run_flash(repo.load_composition("KRSNL_PVTSIM"), 200, 80)
    assert isinstance(result, FlashResult)
    assert result.is_two_phase is False
    assert result.liquid.mole_fraction == pytest.approx(1.0)
    assert result.vapor.mole_fraction == pytest.approx(0.0)


def test_property_rows_keys_present(repo):
    result = flash_service.run_flash(repo.load_composition("KRSNL_PVTSIM"), 50, 20)
    for key, _label in flash_service.PHASE_PROPERTY_ROWS:
        assert key in result.liquid.properties, key


# --- узел «Флэш» в AppState -----------------------------------------------

def test_flash_node_created_on_open(state):
    node = state.flash_node
    assert node is not None
    assert node.kind is NodeKind.FLASH
    assert node.status is NodeStatus.EMPTY
    assert node.upstream == ["composition"]
    assert "P" in node.params and "T" in node.params


def test_flash_status_transitions(state):
    events: list = []
    state.subscribe(lambda: events.append(state.flash_node.status))
    node = state.flash_node

    state.set_flash_running()
    assert node.status is NodeStatus.RUNNING

    fake = flash_service.run_flash(state.active_composition, 200, 80)
    state.set_flash_result(fake)
    assert node.status is NodeStatus.FRESH
    assert node.result is fake

    state.clear_flash()
    assert node.status is NodeStatus.EMPTY
    assert len(events) == 3  # running, result, clear


def test_flash_error(state):
    state.set_flash_error("boom")
    assert state.flash_node.status is NodeStatus.STALE
    assert state.flash_node.error == "boom"
    assert state.flash_node.result is None


def test_composition_change_invalidates_flash(state):
    node = state.flash_node
    state.set_flash_result(flash_service.run_flash(state.active_composition, 200, 80))
    assert node.status is NodeStatus.FRESH
    # любая правка состава помечает посчитанный флэш STALE
    state.normalize_composition()
    assert node.status is NodeStatus.STALE


def test_set_flash_params(state):
    state.set_flash_params(123.0, 45.0)
    assert state.flash_node.params["P"] == 123.0
    assert state.flash_node.params["T"] == 45.0
