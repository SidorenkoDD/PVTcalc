"""
Тесты Фазы 5: сохранение/восстановление сессии — слепок результата флэша
и сериализация `SessionState`. Без DearPyGui, идут в CI.
"""

from pathlib import Path

import pytest

from gui.services import flash_service
from gui.services.model_repository import ModelRepository
from gui.session import SessionState, load_session, save_session

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"


@pytest.fixture
def repo() -> ModelRepository:
    return ModelRepository(db_path=str(MODELS_JSON))


def test_flash_snapshot_roundtrip(repo):
    result = flash_service.run_flash(repo.load_composition("KRSNL_PVTSIM"), 50, 20)
    snap = flash_service.snapshot_flash_result(result)

    # слепок сериализуем в JSON (только числа/строки/списки/словари)
    import json
    json.dumps(snap)

    restored = flash_service.restore_flash_result(snap)
    assert restored.is_two_phase == result.is_two_phase
    assert restored.pressure == pytest.approx(result.pressure)
    assert restored.temperature == pytest.approx(result.temperature)
    assert restored.vapor.mole_fraction == pytest.approx(result.vapor.mole_fraction)
    assert restored.liquid.properties["density"] == pytest.approx(
        float(result.liquid.properties["density"]))


def test_session_v2_workspaces_roundtrip(tmp_path):
    path = str(tmp_path / "sess.json")
    ws = {"open_tabs": ["composition", "flash_1"], "active_tab": "flash_1",
          "flashes": [{"P": 50.0, "T": 20.0, "result": {"is_two_phase": True}}],
          "experiments": []}
    original = SessionState(active_model_id="KRSNL_PVTSIM", window_width=1000,
                            workspaces={"KRSNL_PVTSIM": ws, "OTHER": {"flashes": []}})
    save_session(original, path)
    loaded = load_session(path)

    assert loaded.version == 3
    assert loaded.active_model_id == "KRSNL_PVTSIM"
    assert loaded.workspaces["KRSNL_PVTSIM"]["active_tab"] == "flash_1"
    assert loaded.workspaces["KRSNL_PVTSIM"]["flashes"][0]["P"] == 50.0
    assert "OTHER" in loaded.workspaces  # чужие workspaces не теряются


def test_session_migration_v1_to_v2(tmp_path):
    # файл старого формата: workspace одной модели в плоских полях
    path = tmp_path / "sess.json"
    path.write_text(
        '{"active_model_id": "KRSNL_PVTSIM", "open_tabs": ["composition"],'
        ' "active_tab": "composition",'
        ' "flashes": [{"P": 50.0, "T": 20.0, "result": null}],'
        ' "experiments": [{"kind": "dle", "params": {}, "result": null}]}',
        encoding="utf-8")
    loaded = load_session(str(path))

    assert loaded.version == 3
    ws = loaded.workspaces["KRSNL_PVTSIM"]
    assert ws["open_tabs"] == ["composition"]
    assert ws["flashes"][0]["P"] == 50.0
    assert ws["experiments"][0]["kind"] == "dle"
    # legacy-поля обнулены после миграции
    assert loaded.flashes is None and loaded.open_tabs is None
    assert loaded.experiments is None and loaded.active_tab is None


def test_save_session_stamps_saved_at(tmp_path):
    path = str(tmp_path / "sess.json")
    s = SessionState(active_model_id="X")
    assert s.saved_at is None
    save_session(s, path)
    assert s.saved_at is not None  # проставлено при сохранении
    assert load_session(path).saved_at == s.saved_at


def test_session_missing_file_defaults(tmp_path):
    loaded = load_session(str(tmp_path / "нет.json"))
    assert loaded.active_model_id is None
    assert loaded.open_tabs is None
    assert loaded.flashes is None


def test_session_load_ignores_unknown_keys(tmp_path):
    path = tmp_path / "sess.json"
    path.write_text('{"active_model_id": "X", "bogus_key": 123}', encoding="utf-8")
    loaded = load_session(str(path))
    assert loaded.active_model_id == "X"
    assert not hasattr(loaded, "bogus_key")


def test_session_invalid_shapes_are_normalized(tmp_path):
    path = tmp_path / "session.json"
    path.write_text(
        '{"active_model_id": 42, "window_width": -10, '
        '"window_height": 50000, "workspaces": ["bad"]}',
        encoding="utf-8",
    )

    loaded = load_session(str(path))

    assert loaded.active_model_id is None
    assert loaded.window_width == 320
    assert loaded.window_height == 10000
    assert loaded.workspaces is None


def test_session_non_object_root_falls_back(tmp_path):
    path = tmp_path / "session.json"
    path.write_text("[1, 2, 3]", encoding="utf-8")

    loaded = load_session(str(path))

    assert loaded == SessionState()
