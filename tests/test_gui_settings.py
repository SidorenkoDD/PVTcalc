"""Тесты read-only settings_service: значения совпадают с движком."""

from calc_core.Utils import Constants
from gui.services import settings_service as ss


def test_engine_defaults_match_constants():
    d = ss.engine_defaults()
    assert d["CONSTANT_R"] == Constants.CONSTANT_R
    assert d["TOL_SAT_PRESSURE"] == Constants.TOL_SAT_PRESSURE
    # стандартные условия читаются из StandardConditions
    assert d["STD_P"] == 1.01325
    assert d["STD_T"] == 293.15


def test_schema_and_engine_values_have_same_keys():
    assert set(ss.ALL_KEYS) == set(ss.engine_defaults())
