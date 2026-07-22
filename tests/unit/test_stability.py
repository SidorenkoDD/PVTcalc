"""Тесты теста двухфазной стабильности (`calc_core/PhaseStability/TwoPhaseStabilityTest.py`).

Смесь инвариантных проверок (высокое давление -> стабильно; ниже насыщения ->
нестабильно; детерминизм; отсутствие мутации состава) и характеризационного пина
`S_v`/`S_l` на реальном составе `KRSNL_PVTSIM`. Температура передаётся в Кельвинах
напрямую, без промежуточного преобразования из °C через `Conditions`.

Как и в flash-тесте, `composition.T = t_K` выставляется до конструирования —
иначе температурно-зависимые BRS-параметры a(T) не согласованы с t теста.
"""

import pytest

from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest

REL = 1e-6

T_RES_K = 110 + 273.15  # 383.15 K — пластовая температура KRSNL

# Двухфазная точка (из flash-теста `test_krsnl_two_phase_below_p_sat`): P=100 бар.
P_TWO_PHASE_BAR = 100.0
# Давление заведомо выше насыщения при этой T — однофазная (стабильная) область.
P_STABLE_BAR = 400.0


def _stability_at(composition, p_bar, t_k):
    """Готовит состав к температуре и строит тест стабильности."""
    composition.T = t_k
    test = TwoPhaseStabilityTest(composition, p_bar, t_k)
    test.calculate_phase_stability()
    return test


def test_high_pressure_is_stable(krsnl_composition):
    """Заведомо выше насыщения: система стабильна, K-значения для флэша не выдаются."""
    test = _stability_at(krsnl_composition, P_STABLE_BAR, T_RES_K)

    assert test.stable is True
    assert test.k_values_for_flash is None


def test_below_saturation_is_unstable(krsnl_composition):
    """Ниже насыщения: система нестабильна, выдаётся dict K-значений на все компоненты."""
    test = _stability_at(krsnl_composition, P_TWO_PHASE_BAR, T_RES_K)

    assert test.stable is False
    assert isinstance(test.k_values_for_flash, dict)
    assert len(test.k_values_for_flash) == len(krsnl_composition.composition)


def test_stability_deterministic(krsnl_composition):
    """Повторный вызов `calculate_phase_stability()` даёт тот же вердикт и те же S_v/S_l."""
    krsnl_composition.T = T_RES_K
    test = TwoPhaseStabilityTest(krsnl_composition, P_TWO_PHASE_BAR, T_RES_K)

    test.calculate_phase_stability()
    stable1, s_v1, s_l1 = test.stable, test.S_v, test.S_l
    test.calculate_phase_stability()

    assert test.stable == stable1
    assert test.S_v == s_v1
    assert test.S_l == s_l1


def test_stability_does_not_mutate_composition(krsnl_composition):
    """Тест стабильности не меняет переданный `Composition` (T и мольные доли)."""
    krsnl_composition.T = T_RES_K
    t_before = krsnl_composition.T
    composition_before = dict(krsnl_composition.composition)

    TwoPhaseStabilityTest(krsnl_composition, P_TWO_PHASE_BAR, T_RES_K).calculate_phase_stability()

    assert krsnl_composition.T == t_before
    assert dict(krsnl_composition.composition) == composition_before


def test_stability_characterization(krsnl_composition):
    """Характеризация: пин трейс-сумм Михельсена S_v/S_l при (100 бар, 383.15 K)."""
    test = _stability_at(krsnl_composition, P_TWO_PHASE_BAR, T_RES_K)

    assert test.S_v == pytest.approx(1.3594802695905899, rel=REL)
    assert test.S_l == pytest.approx(0.9999995334454227, rel=REL)
