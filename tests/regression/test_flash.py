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
    assert result.phase_type == "vapor"
    assert result.vapor.mole_fraction == 0.0
    assert result.liquid.mole_fraction == 1.0

    props = result.liquid.properties
    assert props["molecular_ weight"] == pytest.approx(90.7050651693483, rel=REL)
    assert props["molar_volume"] == pytest.approx(232.1824045791991, rel=REL)
    assert props["density"] == pytest.approx(0.39066295886520613, rel=REL)
    assert props["z"] == pytest.approx(1.052745404441804, rel=REL)
    assert props["viscosity"] == pytest.approx(0.05486606329262155, rel=REL)


def test_krsnl_two_phase_below_p_sat(krsnl_composition):
    """KRSNL_PVTSIM, P=100 бар, T=110°C — двухфазная область."""
    conditions = Conditions(100, 110)
    krsnl_composition.T = conditions.t

    result = Flash(krsnl_composition, conditions).calculate()

    assert result.is_two_phase is True
    assert result.vapor.mole_fraction == pytest.approx(0.3233814250596266, rel=REL)
    assert result.liquid.mole_fraction == pytest.approx(0.6766185749403735, rel=REL)

    liquid = result.liquid.properties
    assert liquid["molecular_ weight"] == pytest.approx(122.86557171002929, rel=REL)
    assert liquid["molar_volume"] == pytest.approx(180.39474497264422, rel=REL)
    assert liquid["density"] == pytest.approx(0.6810928540554833, rel=REL)
    assert liquid["z"] == pytest.approx(0.5662664965468325, rel=REL)
    assert liquid["viscosity"] == pytest.approx(0.29192545629531313, rel=REL)

    vapor = result.vapor.properties
    assert vapor["molecular_ weight"] == pytest.approx(23.414879590804794, rel=REL)
    assert vapor["molar_volume"] == pytest.approx(271.03681369721, rel=REL)
    assert vapor["density"] == pytest.approx(0.08639003414850809, rel=REL)
    assert vapor["z"] == pytest.approx(0.8507956645345173, rel=REL)
    assert vapor["viscosity"] == pytest.approx(0.015580193236990185, rel=REL)


def test_przlm_single_phase(przlm_composition):
    """PRRZLM_MDT_TEST, P=200 бар, T=100°C — второй состав, другая молекулярная масса/плотность.

    Эталон обновлён 2026-07-19 по решению автора. Причина — не изменение кода:
    запись PRRZLM_MDT_TEST в models.json была пересохранена через GUI
    (updated_at 2026-07-19T14:56:28), и свойства C7+ пересчитались по
    записанному набору корреляций (pedersen Tc / riazi daubert Pc). Изменились
    critical_temperature/critical_pressure у 24 компонент (например,
    C30 Tc: 872.53 -> 842.92 K) и, следом, коэффициенты EOS a/b/c/d и
    peneloux_correction у всех 35. Автор подтвердил, что верны новые значения:
    это честный пересчёт текущим кодом, тогда как старые лежали в базе с июня.

    Проверено, что дело в данных, а не в коде: на версии models.json из HEAD
    тест проходил со старым эталоном. `molecular_ weight` не изменился —
    сам состав тот же, поехали только EOS-зависимые свойства.

    После итерации R0.1 вход зафиксирован в `tests/fixtures/models.json`:
    дальнейшая работа в GUI с корневой базой этот baseline не сдвигает.
    """
    conditions = Conditions(200, 100)
    przlm_composition.T = conditions.t

    result = Flash(przlm_composition, conditions).calculate()

    assert result.is_two_phase is False
    assert result.phase_type == "ambiguous"
    props = result.liquid.properties
    assert props["molecular_ weight"] == pytest.approx(162.37014542290848, rel=REL)
    assert props["molar_volume"] == pytest.approx(223.02467067381173, rel=REL)
    assert props["density"] == pytest.approx(0.7280367007486159, rel=REL)
    assert props["z"] == pytest.approx(1.4376900642941304, rel=REL)
    assert props["viscosity"] == pytest.approx(0.7167218269379607, rel=REL)
