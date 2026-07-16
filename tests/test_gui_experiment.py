"""
Тесты экспериментов (сервис + узлы AppState). Гоняем DLE (быстрый, ~0.3с);
CCE намеренно не в CI — она тяжёлая (несколько секунд) и капризна к Psat.
"""

import json
from pathlib import Path

import pytest

from gui.app_state import AppState, NodeKind, NodeStatus
from gui.services import experiment_service as exp_svc
from gui.services.model_repository import ModelRepository
from gui.session import SessionState, load_session, save_session

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"

_DLE_PARAMS = {"pressures": [400, 300, 200, 100, 50], "T_c": 100.0, "P_res": 400.0}


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

def test_experiment_types_metadata():
    assert set(exp_svc.EXPERIMENT_TYPES) == {"cce", "dle", "separator"}
    assert exp_svc.EXPERIMENT_TYPES["dle"]["plot_y"] == ["Bo", "Rs"]
    assert exp_svc.EXPERIMENT_TYPES["separator"]["needs_stage_temps"] is True


def test_run_dle_returns_json_table(repo):
    result = exp_svc.run_experiment(repo.load_composition("KRSNL_PVTSIM"),
                                    "dle", _DLE_PARAMS)
    assert result["x"] == "pressure"
    assert "Bo" in result["columns"] and "Rs" in result["columns"]
    assert result["plot_y"] == ["Bo", "Rs"]
    assert len(result["rows"]) > 0
    json.dumps(result)  # JSON-совместимо (NaN → None)


def test_series_for_plot_sorted_no_none(repo):
    result = exp_svc.run_experiment(repo.load_composition("KRSNL_PVTSIM"),
                                    "dle", _DLE_PARAMS)
    xs, ys = exp_svc.series_for_plot(result, "Bo")
    assert len(xs) == len(ys) and len(xs) > 0
    assert xs == sorted(xs)
    assert all(v is not None for v in xs + ys)


def test_run_unknown_kind_raises(repo):
    with pytest.raises(ValueError):
        exp_svc.run_experiment(repo.load_composition("KRSNL_PVTSIM"), "bogus", {})


# --- узлы AppState ---------------------------------------------------------

def test_new_experiment_node(state):
    nid = state.new_experiment("dle", {"pressures": [400, 200], "T_c": 90.0,
                                       "P_res": 400.0})
    node = state.node_by_id(nid)
    assert node.kind is NodeKind.EXPERIMENT
    assert node.params["kind"] == "dle"
    assert node.status is NodeStatus.EMPTY
    assert nid in state.active_variant.open_node_ids
    assert node in state.active_variant.experiment_runs()


def test_composition_change_invalidates_experiment(state):
    nid = state.new_experiment("dle", {"pressures": [400, 200], "T_c": 90.0})
    state.set_node_result(nid, {"columns": ["pressure"], "rows": [[1.0]],
                                "x": "pressure", "plot_y": []})
    assert state.node_by_id(nid).status is NodeStatus.FRESH
    state.normalize_composition()
    assert state.node_by_id(nid).status is NodeStatus.STALE


def test_undo_covers_new_experiment(state):
    nid = state.new_experiment("dle", {"pressures": [400, 200], "T_c": 90.0})
    assert nid in state.active_variant.nodes
    state.undo()
    assert nid not in state.active_variant.nodes


def test_session_experiments_roundtrip(tmp_path):
    path = str(tmp_path / "s.json")
    s = SessionState(
        active_model_id="KRSNL_PVTSIM",
        experiments=[{"kind": "dle",
                      "params": {"kind": "dle", "pressures": [400, 200], "T_c": 90.0},
                      "result": {"columns": ["pressure", "Bo"], "rows": [[400, 1.4]],
                                 "x": "pressure", "plot_y": ["Bo"]}}],
    )
    save_session(s, path)
    loaded = load_session(path)
    assert loaded.experiments[0]["kind"] == "dle"
    assert loaded.experiments[0]["result"]["plot_y"] == ["Bo"]
