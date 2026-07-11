"""Regression-тесты на Flash.calculate().

Это НЕ проверка физической правильности расчёта (нет независимого эталона —
см. CLAUDE.md, раздел про models.json и PVTSim). Тесты фиксируют текущий вывод
движка на реальных составах и ловят непреднамеренные изменения результата при
рефакторинге. Если число законно меняется (осознанно поправили формулу/баг) —
эталон в тесте нужно пересчитать и обновить сознательно, а не потому что тест мешает.
"""

import pytest

from calc_core.Utils.Conditions import Conditions
from calc_core.VLE.Flash import Flash

REL = 1e-6


def test_krsnl_single_phase_above_p_sat(krsnl_composition):
    """KRSNL_PVTSIM, P=250 бар, T=390°C — однофазная область (стабильно)."""
    conditions = Conditions(250, 390)
    krsnl_composition.T = conditions.t

    result = Flash(krsnl_composition, conditions).calculate()

    assert result.is_two_phase is False
    assert result.vapor.mole_fraction == 0.0
    assert result.liquid.mole_fraction == 1.0

    props = result.liquid.properties
    assert props["molecular_ weight"] == pytest.approx(90.7050651693483, rel=REL)
    assert props["molar_volume"] == pytest.approx(232.05312309533375, rel=REL)
    assert props["density"] == pytest.approx(0.3908806050935337, rel=REL)
    assert props["z"] == pytest.approx(1.052740128380821, rel=REL)
    assert props["viscosity"] == pytest.approx(0.05490972013659039, rel=REL)


def test_krsnl_two_phase_below_p_sat(krsnl_composition):
    """KRSNL_PVTSIM, P=100 бар, T=110°C — двухфазная область."""
    conditions = Conditions(100, 110)
    krsnl_composition.T = conditions.t

    result = Flash(krsnl_composition, conditions).calculate()

    assert result.is_two_phase is True
    assert result.vapor.mole_fraction == pytest.approx(0.323361534981212, rel=REL)
    assert result.liquid.mole_fraction == pytest.approx(0.676638465018788, rel=REL)

    liquid = result.liquid.properties
    assert liquid["molecular_ weight"] == pytest.approx(122.86293663059008, rel=REL)
    assert liquid["molar_volume"] == pytest.approx(180.29311593492352, rel=REL)
    assert liquid["density"] == pytest.approx(0.6814621622876452, rel=REL)
    assert liquid["z"] == pytest.approx(0.5662661819485649, rel=REL)
    assert liquid["viscosity"] == pytest.approx(0.2930380292976516, rel=REL)

    vapor = result.vapor.properties
    assert vapor["molecular_ weight"] == pytest.approx(23.414276284794113, rel=REL)
    assert vapor["molar_volume"] == pytest.approx(270.88208995467915, rel=REL)
    assert vapor["density"] == pytest.approx(0.086437151635649, rel=REL)
    assert vapor["z"] == pytest.approx(0.8507888170963235, rel=REL)
    assert vapor["viscosity"] == pytest.approx(0.015582436301816724, rel=REL)


def test_przlm_single_phase(przlm_composition):
    """PRRZLM_MDT_TEST, P=200 бар, T=100°C — второй состав, другая молекулярная масса/плотность."""
    conditions = Conditions(200, 100)
    przlm_composition.T = conditions.t

    result = Flash(przlm_composition, conditions).calculate()

    assert result.is_two_phase is False
    props = result.liquid.properties
    assert props["molecular_ weight"] == pytest.approx(162.37014542290848, rel=REL)
    assert props["molar_volume"] == pytest.approx(233.68588095320453, rel=REL)
    assert props["density"] == pytest.approx(0.694822232137435, rel=REL)
    assert props["z"] == pytest.approx(1.5072650822412388, rel=REL)
    assert props["viscosity"] == pytest.approx(0.4672075164503361, rel=REL)
