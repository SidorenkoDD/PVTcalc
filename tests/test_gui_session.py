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


def test_session_save_load_new_fields(tmp_path):
    path = str(tmp_path / "sess.json")
    original = SessionState(
        active_model_id="KRSNL_PVTSIM",
        window_width=1000,
        window_height=700,
        composition_window_open=True,
        flash={"P": 50.0, "T": 20.0, "result": {"is_two_phase": True}},
    )
    save_session(original, path)
    loaded = load_session(path)

    assert loaded.active_model_id == "KRSNL_PVTSIM"
    assert loaded.window_width == 1000
    assert loaded.composition_window_open is True
    assert loaded.flash["P"] == 50.0
    assert loaded.flash["result"]["is_two_phase"] is True


def test_session_missing_file_defaults(tmp_path):
    loaded = load_session(str(tmp_path / "нет.json"))
    assert loaded.active_model_id is None
    assert loaded.composition_window_open is False
    assert loaded.flash is None


def test_session_load_ignores_unknown_keys(tmp_path):
    path = tmp_path / "sess.json"
    path.write_text('{"active_model_id": "X", "bogus_key": 123}', encoding="utf-8")
    loaded = load_session(str(path))
    assert loaded.active_model_id == "X"
    assert not hasattr(loaded, "bogus_key")
