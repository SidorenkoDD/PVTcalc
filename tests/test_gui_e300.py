"""
Тесты экспорта/импорта Eclipse 300 (.inc): нативные единицы METRIC и
round-trip сохранения свойств. Без DearPyGui (CI-friendly).
"""

import json
from pathlib import Path

import pytest

from calc_core.Composition.Composition import Composition
from calc_core.Utils.E300Export import E300Exporter
from calc_core.Utils.E300Import import parse_e300
from gui.services import export_service as exp_out_svc
from gui.services import project_service as proj_svc

MODELS_JSON = Path(__file__).resolve().parents[1] / "models.json"


@pytest.fixture
def comp() -> Composition:
    db = Composition.from_db(str(MODELS_JSON))
    return getattr(db, "KRSNL_PVTSIM")


# --- экспорт ---------------------------------------------------------------

def test_e300_export_native_units(comp):
    s = E300Exporter(comp, "MPR").export_string()
    # структура блоков
    for kw in ("EOS", "NCOMPS", "CNAMES", "MW", "TCRIT", "PCRIT", "ZI", "BIC"):
        assert f"\n{kw}\n" in s
    assert "\nMPR\n" in s
    # нативные единицы: MW=г/моль (16.043, не 0.016), Pc=бар (46.0, не 460)
    assert "16.043" in s
    assert "\n46.0\n" in s


def test_e300_export_eos_keyword_validated(comp):
    with pytest.raises(ValueError):
        E300Exporter(comp, "PR78")  # не MPR/SRK


def test_export_service_writes_file(comp, tmp_path):
    out = str(tmp_path / "deck")
    written = exp_out_svc.export_model(comp, "e300", out, eos_keyword="SRK")
    assert written.endswith(".inc")
    text = Path(written).read_text(encoding="utf-8")
    assert "\nSRK\n" in text


# --- парсер ----------------------------------------------------------------

def test_parse_e300_roundtrip(comp, tmp_path):
    deck = str(tmp_path / "d.inc")
    E300Exporter(comp, "MPR").export(deck)
    parsed = parse_e300(deck)
    names = list(comp.composition.keys())
    assert parsed["names"] == names
    assert parsed["eos"] == "MPR"
    # zi и свойства совпали с исходными (нативные единицы, без конверсий)
    assert parsed["zi"]["C1"] == pytest.approx(comp.composition["C1"])
    assert parsed["props"]["critical_pressure"]["C1"] == pytest.approx(46.0)
    assert parsed["props"]["molar_mass"]["C1"] == pytest.approx(16.043)


def test_parse_e300_requires_cnames_zi(tmp_path):
    bad = tmp_path / "bad.inc"
    bad.write_text("EOS\n--\nMPR\n/\n", encoding="utf-8")
    with pytest.raises(ValueError):
        parse_e300(str(bad))


# --- импорт (round-trip во временную БД) -----------------------------------

def test_import_e300_preserves_properties(comp, tmp_path):
    deck = str(tmp_path / "rt.inc")
    E300Exporter(comp, "SRK").export(deck)
    db_path = str(tmp_path / "models.json")

    res = proj_svc.import_e300(db_path, deck)
    assert res["recognized"] and not res["unrecognized"]
    assert res["deck_eos"] == "SRK"

    saved = json.loads(Path(db_path).read_text(encoding="utf-8"))[res["model_id"]]
    cd = saved["composition_data"]
    # свойства дека перенесены как есть (нативные единицы сохранены)
    assert cd["critical_pressure"]["C1"] == pytest.approx(
        comp.composition_data["critical_pressure"]["C1"])
    assert cd["critical_temperature"]["C1"] == pytest.approx(
        comp.composition_data["critical_temperature"]["C1"])
    assert cd["molar_mass"]["C1"] == pytest.approx(
        comp.composition_data["molar_mass"]["C1"])
    # модель сохранена на живом EOS Брусиловского независимо от ключа дека
    assert saved["eos"] == "BRSEOS"
