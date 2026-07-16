"""
Тесты фреймворк-независимого слоя GUI (`gui.app_state` + `gui.services`).

Не импортируют DearPyGui — проверяют только состояние и загрузку моделей
из `models.json`, поэтому идут в CI без опциональной группы `[gui]`.
"""

from pathlib import Path

import pytest

from gui.app_state import AppState, NodeKind, NodeStatus
from gui.services.model_repository import ModelRepository

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"


@pytest.fixture
def repo() -> ModelRepository:
    return ModelRepository(db_path=str(MODELS_JSON))


def test_repository_lists_models(repo):
    summaries = repo.list_models()
    ids = {s.model_id for s in summaries}
    assert {"KRSNL_PVTSIM", "PRRZLM_MDT_TEST"} <= ids
    krsnl = next(s for s in summaries if s.model_id == "KRSNL_PVTSIM")
    assert krsnl.n_components == 40
    assert krsnl.eos == "PREOS"


def test_repository_missing_file_returns_empty(tmp_path):
    repo = ModelRepository(db_path=str(tmp_path / "нет.json"))
    assert repo.list_models() == []


def test_repository_load_composition(repo):
    comp = repo.load_composition("KRSNL_PVTSIM")
    assert len(comp.composition) == 40
    # сохранённые свойства инжектированы без пересчёта корреляций
    assert comp.composition_data["critical_temperature"]["C1"] == pytest.approx(190.6)


def test_repository_load_unknown_raises(repo):
    with pytest.raises(KeyError):
        repo.load_composition("НЕТ_ТАКОЙ")


def test_app_state_refresh_and_open(repo):
    state = AppState(repo)
    notifications: list = []
    state.subscribe(lambda: notifications.append(state.active_model_id))

    state.refresh_model_list()
    assert "KRSNL_PVTSIM" in state.models
    assert state.models["KRSNL_PVTSIM"].loaded is False

    state.open_model("KRSNL_PVTSIM")
    assert state.active_model_id == "KRSNL_PVTSIM"
    assert state.active_variant_id == "base"
    assert state.active_composition is not None
    assert len(state.active_composition.composition) == 40

    node = state.active_variant.nodes["composition"]
    assert node.kind is NodeKind.COMPOSITION
    assert node.status is NodeStatus.FRESH

    # refresh + open породили хотя бы два уведомления наблюдателя
    assert len(notifications) >= 2


def test_app_state_reopen_reuses_loaded_variant(repo):
    state = AppState(repo)
    state.refresh_model_list()
    state.open_model("KRSNL_PVTSIM")
    first_comp = state.active_composition

    # повторный refresh не должен сбрасывать уже загруженный вариант
    state.refresh_model_list()
    assert state.models["KRSNL_PVTSIM"].loaded is True
    state.open_model("KRSNL_PVTSIM")
    assert state.active_composition is first_comp


def test_open_unknown_model_raises(repo):
    state = AppState(repo)
    state.refresh_model_list()
    with pytest.raises(KeyError):
        state.open_model("НЕТ")
