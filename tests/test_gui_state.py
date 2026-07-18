"""
Тесты фреймворк-независимого слоя GUI (`gui.app_state` + `gui.services`).

Не импортируют DearPyGui — проверяют только состояние и загрузку моделей
из `models.json`, поэтому идут в CI без опциональной группы `[gui]`.
"""

from pathlib import Path

import pytest

from gui.app_state import (
    AppState,
    NodeKind,
    NodeStatus,
    StateChangeKind,
)
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


def test_app_state_exposes_corrupt_store_without_crashing(tmp_path):
    db_path = tmp_path / "models.json"
    db_path.write_text('{"broken":', encoding="utf-8")
    state = AppState(ModelRepository(db_path=str(db_path)))
    notifications: list[None] = []
    state.subscribe(lambda: notifications.append(None))

    state.refresh_model_list()

    assert state.models == {}
    assert state.model_list_error is not None
    assert "строка 1, столбец 11" in state.model_list_error
    assert notifications == [None]

    db_path.write_text("{}", encoding="utf-8")
    state.refresh_model_list()
    assert state.model_list_error is None


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


def test_state_change_events_describe_render_scope(repo):
    state = AppState(repo)
    changes = []
    state.subscribe_changes(changes.append)

    state.refresh_model_list()
    state.open_model("KRSNL_PVTSIM")
    state.open_node("composition")
    flash_id = state.new_flash_run(100.0, 20.0)
    state.set_node_running(flash_id)

    assert [change.kind for change in changes] == [
        StateChangeKind.MODEL_LIST,
        StateChangeKind.WORKSPACE,
        StateChangeKind.WORKSPACE,
        StateChangeKind.WORKSPACE,
        StateChangeKind.NODE,
    ]
    assert changes[-1].node_ref is not None
    assert changes[-1].node_ref.node_id == flash_id


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


# --- навигация между экранами (Projects / New fluid / Workspace) ------------

def test_starts_on_projects_screen(repo):
    state = AppState(repo)
    assert state.current_screen == "projects"


def test_enter_model_switches_to_workspace(repo):
    state = AppState(repo)
    state.refresh_model_list()
    events: list = []
    state.subscribe(lambda: events.append(state.current_screen))
    state.enter_model("KRSNL_PVTSIM")
    assert state.current_screen == "workspace"
    assert state.active_model_id == "KRSNL_PVTSIM"
    assert state.models["KRSNL_PVTSIM"].loaded is True
    assert events[-1] == "workspace"


def test_show_projects_keeps_active_model(repo):
    state = AppState(repo)
    state.refresh_model_list()
    state.enter_model("KRSNL_PVTSIM")
    state.show_projects()
    assert state.current_screen == "projects"
    # активная модель не сбрасывается — «Continue last» и workspace живы
    assert state.active_model_id == "KRSNL_PVTSIM"


def test_enter_unknown_model_keeps_screen(repo):
    state = AppState(repo)
    state.refresh_model_list()
    with pytest.raises(KeyError):
        state.enter_model("НЕТ")
    assert state.current_screen == "projects"


def test_model_summary_attached(repo):
    state = AppState(repo)
    state.refresh_model_list()
    s = state.models["KRSNL_PVTSIM"].summary
    assert s is not None
    assert s.results_brief and s.results_brief[0]["module"] == "Flash"


def test_node_ref_resolves_after_switching_models(repo):
    state = AppState(repo)
    state.refresh_model_list()
    state.open_model("KRSNL_PVTSIM")
    nid = state.new_flash_run(50, 20)
    ref = state.node_ref(nid)
    state.open_model("PRRZLM_MDT_TEST")
    marker = {"done": True}
    state.set_node_result(ref, marker)
    assert state.node_by_ref(ref).result is marker
    assert state.node_by_id(nid) is None
