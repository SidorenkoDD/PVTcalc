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

MODELS_JSON = Path(__file__).resolve().parent / "fixtures" / "models.json"


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
    assert result.quality_status == "ok"
    assert result.diagnostics.warnings == ()


def test_run_flash_single_phase(repo):
    result = flash_service.run_flash(repo.load_composition("KRSNL_PVTSIM"), 200, 80)
    assert isinstance(result, FlashResult)
    assert result.is_two_phase is False
    assert result.liquid.mole_fraction == pytest.approx(1.0)
    assert result.vapor.mole_fraction == pytest.approx(0.0)
    assert result.quality_status == "ok"


def test_property_rows_keys_present(repo):
    result = flash_service.run_flash(repo.load_composition("KRSNL_PVTSIM"), 50, 20)
    for key, _label in flash_service.PHASE_PROPERTY_ROWS:
        assert key in result.liquid.properties, key


# --- узлы «Флэш» (история) в AppState -------------------------------------

def test_new_flash_run_creates_node_and_tab(state):
    nid = state.new_flash_run(50, 20)
    node = state.node_by_id(nid)
    assert node.kind is NodeKind.FLASH
    assert node.status is NodeStatus.EMPTY
    assert node.upstream == ["composition"]
    assert node.params["P"] == 50 and node.params["T"] == 20
    # открылся вкладкой и стал активным
    assert nid in state.active_variant.open_node_ids
    assert state.active_variant.active_node_id == nid


def test_history_multiple_runs_ordered(state):
    a = state.new_flash_run(50, 20)
    b = state.new_flash_run(200, 80)
    assert [r.node_id for r in state.active_variant.flash_runs()] == [a, b]


def test_flash_status_transitions(state):
    nid = state.new_flash_run(200, 80)
    events: list = []
    state.subscribe(lambda: events.append(1))

    state.set_flash_running(nid)
    assert state.node_by_id(nid).status is NodeStatus.RUNNING

    fake = flash_service.run_flash(state.active_composition, 200, 80)
    state.set_flash_result(nid, fake)
    assert state.node_by_id(nid).status is NodeStatus.FRESH
    assert state.node_by_id(nid).result is fake

    state.reset_flash(nid)
    assert state.node_by_id(nid).status is NodeStatus.EMPTY
    assert len(events) == 3  # running, result, reset


def test_flash_error(state):
    nid = state.new_flash_run(200, 80)
    state.set_flash_error(nid, "boom")
    node = state.node_by_id(nid)
    assert node.status is NodeStatus.STALE
    assert node.error == "boom"
    assert node.result is None


def test_composition_change_invalidates_all_flashes(state):
    a = state.new_flash_run(50, 20)
    b = state.new_flash_run(200, 80)
    state.set_flash_result(a, flash_service.run_flash(state.active_composition, 50, 20))
    state.set_flash_result(b, flash_service.run_flash(state.active_composition, 200, 80))
    assert state.node_by_id(a).status is NodeStatus.FRESH
    # любая правка состава помечает ВСЕ посчитанные флэши STALE
    state.normalize_composition()
    assert state.node_by_id(a).status is NodeStatus.STALE
    assert state.node_by_id(b).status is NodeStatus.STALE


def test_set_flash_params(state):
    nid = state.new_flash_run()
    state.set_flash_params(nid, 123.0, 45.0)
    assert state.node_by_id(nid).params["P"] == 123.0
    assert state.node_by_id(nid).params["T"] == 45.0


def test_close_node_updates_tabs(state):
    a = state.new_flash_run(50, 20)
    b = state.new_flash_run(200, 80)
    assert state.active_variant.active_node_id == b
    state.close_node(b)
    assert b not in state.active_variant.open_node_ids
    assert state.active_variant.active_node_id == a  # активной стала соседняя


def test_duplicate_flash(state):
    a = state.new_flash_run(50, 20)
    b = state.duplicate_flash(a)
    assert b != a
    node = state.node_by_id(b)
    assert node.params["P"] == 50 and node.params["T"] == 20
    assert node.status is NodeStatus.EMPTY  # копия не запущена
    assert b in state.active_variant.open_node_ids


def test_rename_node_sets_and_clears_label(state):
    a = state.new_flash_run(50, 20)
    state.rename_node(a, "  bubble @ res  ")
    assert state.node_by_id(a).params["label"] == "bubble @ res"
    state.rename_node(a, "   ")  # пусто -> сброс
    assert "label" not in state.node_by_id(a).params


def test_delete_node_removes_from_history_and_tabs(state):
    a = state.new_flash_run(50, 20)
    b = state.new_flash_run(200, 80)
    state.delete_node(a)
    assert a not in state.active_variant.nodes
    assert a not in state.active_variant.open_node_ids
    assert b in state.active_variant.nodes


def test_delete_composition_is_blocked(state):
    state.open_node("composition")
    state.delete_node("composition")
    assert "composition" in state.active_variant.nodes  # состав не удаляется


def test_open_compare_orders_by_creation(state):
    a = state.new_flash_run(50, 20)
    b = state.new_flash_run(200, 80)
    nid = state.open_compare([b, a])  # выбор в обратном порядке
    node = state.node_by_id(nid)
    assert node.kind is NodeKind.COMPARE
    assert node.params["members"] == [a, b]  # порядок создания, не выбора
    assert nid in state.active_variant.open_node_ids


def test_open_compare_needs_two(state):
    a = state.new_flash_run(50, 20)
    assert state.open_compare([a]) is None  # одного мало
