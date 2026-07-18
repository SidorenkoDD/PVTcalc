"""Тесты settings_service: дефолты движка + round-trip сохранения. Без DPG."""

import json
from pathlib import Path

from calc_core.Utils import Constants
from gui.services import settings_service as ss


def test_engine_defaults_match_constants():
    d = ss.engine_defaults()
    assert d["CONSTANT_R"] == Constants.CONSTANT_R
    assert d["TOL_SAT_PRESSURE"] == Constants.TOL_SAT_PRESSURE
    # стандартные условия читаются из StandardConditions
    assert d["STD_P"] == 1.01325
    assert d["STD_T"] == 293.15


def test_save_load_roundtrip(tmp_path):
    path = str(tmp_path / "gui_settings.json")
    values = ss.engine_defaults()
    values["CONSTANT_R"] = 8.314462618
    ss.save_settings(values, path)
    # на диске только ключи схемы
    saved = json.loads(Path(path).read_text(encoding="utf-8"))
    assert set(saved).issubset(set(ss.ALL_KEYS))
    # загрузка накладывает правки поверх дефолтов
    loaded = ss.load_settings(path)
    assert loaded["CONSTANT_R"] == 8.314462618


def test_load_missing_file_returns_defaults(tmp_path):
    loaded = ss.load_settings(str(tmp_path / "nope.json"))
    assert loaded["CONSTANT_R"] == Constants.CONSTANT_R


def test_load_corrupt_file_falls_back(tmp_path):
    p = tmp_path / "bad.json"
    p.write_text("{ not json", encoding="utf-8")
    loaded = ss.load_settings(str(p))
    assert loaded["CONSTANT_R"] == Constants.CONSTANT_R
