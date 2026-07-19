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


def test_run_dle_returns_rich_json_result(repo):
    result = exp_svc.run_experiment(repo.load_composition("KRSNL_PVTSIM"),
                                    "dle", _DLE_PARAMS)
    assert result["x"] == "pressure"
    assert "Bo" in result["columns"] and "Rs" in result["columns"]
    assert result["plot_y"] == ["Bo", "Rs"]
    assert result["main_columns"][0] == "pressure"
    assert set(["Bo", "Rs"]).issubset(result["charts"])  # больше одной кривой
    assert "pressure" not in result["plot_all"]  # x не строим как кривую
    assert len(result["rows"]) > 0
    json.dumps(result)  # JSON-совместимо (NaN → None), включая stages


def test_run_dle_has_per_stage_compositions(repo):
    result = exp_svc.run_experiment(repo.load_composition("KRSNL_PVTSIM"),
                                    "dle", _DLE_PARAMS)
    stages = result["stages"]
    assert len(stages) == len(result["rows"])
    first = stages[0]
    assert isinstance(first["liquid"], dict) and "C1" in first["liquid"]
    assert isinstance(first["vapor"], dict)


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


def test_lab_data_rows_are_typed_and_editable(state):
    nid = state.new_experiment("dle", {"pressures": [400, 200]})
    columns = exp_svc.EXPERIMENT_TYPES["dle"]["lab_columns"]

    state.add_lab_data_row(nid, columns)
    state.set_lab_data_value(nid, 0, 0, 250.0)
    state.set_lab_data_value(nid, 0, 1, 1.23)
    state.append_lab_data_rows(nid, columns, [[200.0, 1.1], [150.0]])
    node = state.node_by_id(nid)
    assert node.params["lab_data"] == {
        "schema_version": 1,
        "columns": columns,
        "rows": [[250.0, 1.23, None, None, None],
                 [200.0, 1.1, None, None, None],
                 [150.0, None, None, None, None]],
    }

    state.add_lab_data_row(nid, columns)
    state.remove_lab_data_row(nid)
    assert len(node.params["lab_data"]["rows"]) == 3
    state.clear_lab_data(nid)
    assert node.params["lab_data"]["rows"] == []


def test_paste_lab_data_rows_starts_at_cell_and_extends_table(state):
    nid = state.new_experiment("dle", {"pressures": [400, 200]})
    columns = exp_svc.EXPERIMENT_TYPES["dle"]["lab_columns"]
    state.append_lab_data_rows(nid, columns, [[400.0, 1.4]])

    state.paste_lab_data_rows(
        nid, columns, start_row=0,
        new_rows=[[None, 1.3], [None, 1.2], [None, 1.1]])

    assert state.node_by_id(nid).params["lab_data"]["rows"] == [
        [400.0, 1.3, None, None, None],
        [None, 1.2, None, None, None],
        [None, 1.1, None, None, None],
    ]


def test_lab_data_edits_support_undo_and_redo(state):
    nid = state.new_experiment("dle", {"pressures": [400, 200]})
    columns = exp_svc.EXPERIMENT_TYPES["dle"]["lab_columns"]

    state.add_lab_data_row(nid, columns)
    state.set_lab_data_value(nid, 0, 0, 350.0)
    assert state.node_by_id(nid).params["lab_data"]["rows"][0][0] == 350.0

    state.undo()
    assert state.node_by_id(nid).params["lab_data"]["rows"][0][0] is None
    state.redo()
    assert state.node_by_id(nid).params["lab_data"]["rows"][0][0] == 350.0

    state.paste_lab_data_rows(nid, columns, 0,
                              [[400.0, 1.4], [300.0, 1.2]])
    state.undo()
    assert state.node_by_id(nid).params["lab_data"]["rows"] == [
        [350.0, None, None, None, None],
    ]
    state.redo()
    assert state.node_by_id(nid).params["lab_data"]["rows"][:2] == [
        [400.0, 1.4, None, None, None],
        [300.0, 1.2, None, None, None],
    ]


def test_extra_experiment_chart_is_stateful_and_undoable(state):
    nid = state.new_experiment("dle", {"pressures": [400, 200]})
    node = state.node_by_id(nid)
    node.result = {"charts": ["Bo"], "plot_all": ["Bo", "Rs"]}

    assert state.add_experiment_chart(nid, "Bo") is False
    assert state.add_experiment_chart(nid, "Rs") is True
    assert state.node_by_id(nid).params["extra_charts"] == ["Rs"]
    state.undo()
    assert "extra_charts" not in state.node_by_id(nid).params
    state.redo()
    assert state.node_by_id(nid).params["extra_charts"] == ["Rs"]


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
    exp = {"kind": "dle",
           "params": {"kind": "dle", "pressures": [400, 200], "T_c": 90.0,
                      "lab_data": {"schema_version": 1,
                                   "columns": ["pressure", "Bo"],
                                   "rows": [[400.0, 1.4]]}},
           "result": {"columns": ["pressure", "Bo"], "rows": [[400, 1.4]],
                      "x": "pressure", "plot_y": ["Bo"]}}
    s = SessionState(active_model_id="KRSNL_PVTSIM",
                     workspaces={"KRSNL_PVTSIM": {"experiments": [exp]}})
    save_session(s, path)
    loaded = load_session(path)
    ws = loaded.workspaces["KRSNL_PVTSIM"]
    assert ws["experiments"][0]["kind"] == "dle"
    assert ws["experiments"][0]["result"]["plot_y"] == ["Bo"]
    assert ws["experiments"][0]["params"]["lab_data"]["rows"] == [[400.0, 1.4]]
