"""
Тесты фазовой огибающей (сервис + узлы AppState).

Движок в CI гоняем на узком/грубом диапазоне (~1 c). Критическая точка (и
особые точки) намеренно НЕ считаются — их расчёт в движке сейчас некорректен
(правка автора 2026-07-17), вернёмся к ним отдельно.
"""

import json
from pathlib import Path

import pytest

from gui.app_state import AppState, NodeKind, NodeStatus
from gui.services import phase_envelope_service as pe_svc
from gui.services.model_repository import ModelRepository
from gui.session import SessionState, load_session, save_session

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"

# узкий/грубый диапазон — быстрый прогон структуры
_NARROW = {"t_min_c": 90.0, "t_max_c": 110.0, "t_step_c": 20.0,
           "p_max_bar": 500.0, "T_res_c": 100.0}


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

def test_default_envelope_params_sane():
    comp = ModelRepository(db_path=str(MODELS_JSON)).load_composition("KRSNL_PVTSIM")
    p = pe_svc.default_envelope_params(comp)
    assert p["t_max_c"] > p["t_min_c"]
    assert p["t_step_c"] > 0
    assert p["p_max_bar"] > 0
    # пластовая температура центрирует окно (T_res = comp.T - 273.15)
    assert p["T_res_c"] == pytest.approx(comp.T - 273.15, abs=1e-6)


def test_run_envelope_ssm_structure_json(repo):
    result = pe_svc.run_envelope(repo.load_composition("KRSNL_PVTSIM"), _NARROW)
    assert result["method"] == "ssm"
    # кривые: три ветки, равные длины T/P, без None внутри
    assert set(result["curves"]) == {"bubble", "dew_upper", "dew_lower"}
    for c in result["curves"].values():
        assert len(c["T"]) == len(c["P"])
        assert all(v is not None for v in c["T"] + c["P"])
    # сшитая линия огибающей (одной кривой)
    env = result["envelope"]
    assert len(env["T"]) == len(env["P"]) and len(env["T"]) > 0
    # таблица покрывает всю сетку температур
    assert result["table"]["columns"][0] == "Temp_C"
    assert len(result["table"]["rows"]) >= 2
    # крит. точка НЕ считается (убрана по правке автора)
    assert "critical" not in result
    # пластовый маркер
    assert result["reservoir"]["T_c"] == 100.0
    assert result["reservoir"]["P_sat"] is not None
    json.dumps(result)  # JSON-совместимо (NaN → None)


def test_joined_envelope_upper_then_lower():
    # верхняя ветка по возрастанию T + нижняя (dew_lower) по убыванию T
    temps = [10.0, 20.0, 30.0]
    raw = {"bubble": [100.0, 120.0, None], "dew_upper": [None, None, 140.0],
           "dew_lower": [50.0, 60.0, None]}
    env = pe_svc._joined_envelope(temps, raw)
    assert env["T"] == [10.0, 20.0, 30.0, 20.0, 10.0]
    assert env["P"] == [100.0, 120.0, 140.0, 60.0, 50.0]


def test_run_envelope_has_bubble_branch(repo):
    # KRSNL_PVTSIM при 90-110°C — нефтяная область: bubble-ветка непуста
    result = pe_svc.run_envelope(repo.load_composition("KRSNL_PVTSIM"), _NARROW)
    assert len(result["curves"]["bubble"]["T"]) > 0


def test_run_grid_scan_json(repo):
    # grid — скан стабильности P×T; малая сетка для скорости
    params = pe_svc.default_envelope_params(repo.load_composition("KRSNL_PVTSIM"))
    params.update({"method": "grid", "grid_t_max_c": 200.0, "grid_p_max_bar": 400.0,
                   "grid_t_points": 5, "grid_p_points": 5})
    result = pe_svc.run_envelope(repo.load_composition("KRSNL_PVTSIM"), params)
    assert result["method"] == "grid"
    uns, sta = result["grid"]["unstable"], result["grid"]["stable"]
    # 5x5=25 точек распределены между двумя множествами
    assert len(uns["T"]) + len(sta["T"]) == 25
    assert result["table"]["columns"] == ["Temp_C", "Pressure_bar", "Two_phase"]
    # для grid Psat не интерполируем — маркер без давления
    assert result["reservoir"]["P_sat"] is None
    json.dumps(result)


def test_reservoir_psat_interpolation():
    temps = [90.0, 110.0]
    primary = [180.0, 200.0]
    # ровно посередине — линейная интерполяция
    assert pe_svc._interp_reservoir_psat(temps, primary, 100.0) == pytest.approx(190.0)
    # на узле сетки — точное значение
    assert pe_svc._interp_reservoir_psat(temps, primary, 90.0) == pytest.approx(180.0)
    # вне диапазона — не экстраполируем
    assert pe_svc._interp_reservoir_psat(temps, primary, 50.0) is None
    assert pe_svc._interp_reservoir_psat(temps, primary, 150.0) is None
    # < 2 валидных точек — None
    assert pe_svc._interp_reservoir_psat([90.0], [180.0], 90.0) is None
    assert pe_svc._interp_reservoir_psat(temps, [180.0, None], 100.0) is None


def test_service_does_not_import_new_methodv2():
    # расчёт идёт только через SSM/фасад — new_methodv2 не импортируется и не
    # вызывается (см. правку автора; упоминание в докстринге допустимо)
    src = Path(pe_svc.__file__).read_text(encoding="utf-8")
    import_lines = [ln for ln in src.splitlines()
                    if ln.lstrip().startswith(("import ", "from "))]
    assert not any("new_methodv2" in ln for ln in import_lines)
    assert "SaturationPressure(" not in src  # нет вызова солвера из new_methodv2


# --- узлы AppState ---------------------------------------------------------

def test_new_envelope_node(state):
    defaults = pe_svc.default_envelope_params(state.active_composition)
    nid = state.new_envelope(defaults)
    node = state.node_by_id(nid)
    assert node.kind is NodeKind.PHASE_ENVELOPE
    assert node.status is NodeStatus.EMPTY
    assert node.upstream == ["composition"]
    assert nid in state.active_variant.open_node_ids
    assert node in state.active_variant.envelope_runs()
    assert nid == "env_1"


def test_composition_change_invalidates_envelope(state):
    nid = state.new_envelope(pe_svc.default_envelope_params(state.active_composition))
    state.set_node_result(nid, {"curves": {}, "table": {"columns": [], "rows": []}})
    assert state.node_by_id(nid).status is NodeStatus.FRESH
    state.normalize_composition()
    assert state.node_by_id(nid).status is NodeStatus.STALE


def test_undo_covers_new_envelope(state):
    nid = state.new_envelope(pe_svc.default_envelope_params(state.active_composition))
    assert nid in state.active_variant.nodes
    state.undo()
    assert nid not in state.active_variant.nodes


def test_delete_envelope_node(state):
    nid = state.new_envelope(pe_svc.default_envelope_params(state.active_composition))
    state.delete_node(nid)
    assert nid not in state.active_variant.nodes
    assert nid not in state.active_variant.open_node_ids


def test_session_envelopes_roundtrip(tmp_path):
    path = str(tmp_path / "s.json")
    env = {"params": {"t_min_c": -50.0, "t_max_c": 150.0, "t_step_c": 10.0,
                      "p_max_bar": 700.0, "T_res_c": 100.0},
           "result": {"curves": {"bubble": {"T": [90.0], "P": [180.0],
                                            "label": "Bubble point"}},
                      "reservoir": {"T_c": 100.0, "P_sat": 178.0},
                      "table": {"columns": ["Temp_C"], "rows": [[90.0]]}}}
    s = SessionState(active_model_id="KRSNL_PVTSIM",
                     workspaces={"KRSNL_PVTSIM": {"envelopes": [env]}})
    save_session(s, path)
    loaded = load_session(path)
    ws = loaded.workspaces["KRSNL_PVTSIM"]
    assert ws["envelopes"][0]["params"]["p_max_bar"] == 700.0
    assert ws["envelopes"][0]["result"]["reservoir"]["P_sat"] == 178.0
