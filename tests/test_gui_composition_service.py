"""
Тесты Фазы 2: сервис операций над составом (`gui.services.composition_service`)
и команды корневого узла в `AppState`. Без DearPyGui — идут в CI.
"""

from pathlib import Path

import pytest

from gui.app_state import AppState, NodeStatus
from gui.services import composition_service as svc
from gui.services.model_repository import ModelRepository

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"


@pytest.fixture
def composition():
    repo = ModelRepository(db_path=str(MODELS_JSON))
    return repo.load_composition("KRSNL_PVTSIM")


@pytest.fixture
def state():
    repo = ModelRepository(db_path=str(MODELS_JSON))
    st = AppState(repo)
    st.refresh_model_list()
    st.open_model("KRSNL_PVTSIM")
    return st


# --- сервис ----------------------------------------------------------------

def test_options_nonempty_and_defaults_valid():
    assert set(svc.EOS_OPTIONS) == {"PREOS", "SRKEOS", "BRSEOS"}
    for prop, method in svc.DEFAULT_C7_CORRELATIONS.items():
        assert method in svc.CORRELATION_OPTIONS[prop], (prop, method)


def test_has_c7_plus(composition):
    assert svc.has_c7_plus(composition) is True


def test_set_eos(composition):
    svc.set_eos(composition, "BRSEOS")
    assert composition.eos_name.value == "BRSEOS"


def test_recalculate_applies_correlation(composition):
    c7 = next(n for n, f in composition.composition_data["c7_plus_flag"].items() if f)
    corr = svc.default_correlations()
    corr["critical_temperature"] = "pedersen"
    svc.recalculate(composition, corr)
    tc_pedersen = composition.composition_data["critical_temperature"][c7]
    corr["critical_temperature"] = "twu"
    svc.recalculate(composition, corr)
    tc_twu = composition.composition_data["critical_temperature"][c7]
    assert tc_pedersen != tc_twu


def test_edit_zi_and_normalize(composition):
    svc.edit_zi(composition, "C1", 0.5)
    assert composition.composition["C1"] == 0.5
    svc.normalize(composition)
    assert svc.sum_zi(composition) == pytest.approx(1.0)


def test_edit_property(composition):
    svc.edit_property(composition, "C1", "critical_temperature", 191.0)
    assert composition.composition_data["critical_temperature"]["C1"] == 191.0


def test_edit_and_get_bip(composition):
    svc.edit_bip(composition, "N2", "C1", 0.05)
    assert svc.get_bip(composition, "N2", "C1") == pytest.approx(0.05)
    # симметрия
    assert svc.get_bip(composition, "C1", "N2") == pytest.approx(0.05)


# --- команды AppState ------------------------------------------------------

def test_set_correlation_marks_stale_then_recalc_fresh(state):
    node = state.composition_node
    assert node.status is NodeStatus.FRESH
    state.set_correlation("critical_temperature", "twu")
    assert node.status is NodeStatus.STALE
    assert node.params["correlations"]["critical_temperature"] == "twu"
    state.recalculate_composition()
    assert node.status is NodeStatus.FRESH
    assert node.error is None


def test_set_eos_updates_node_and_composition(state):
    state.set_composition_eos("SRKEOS")
    assert state.active_composition.eos_name.value == "SRKEOS"
    assert state.composition_node.params["eos"] == "SRKEOS"


def test_immediate_edits_do_not_notify(state):
    calls: list = []
    state.subscribe(lambda: calls.append(1))
    state.edit_zi("C1", 0.4)
    state.edit_component_property("C1", "molar_mass", 16.1)
    state.edit_bip("N2", "C1", 0.02)
    assert calls == []  # immediate-правки не дёргают наблюдателя


def test_normalize_notifies(state):
    calls: list = []
    state.subscribe(lambda: calls.append(1))
    state.edit_zi("C1", 0.4)
    state.normalize_composition()
    assert len(calls) == 1
    assert sum(state.active_composition.composition.values()) == pytest.approx(1.0)
