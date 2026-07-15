"""Тесты EOS Брусиловского (`calc_core/EOS/BrusilovskiyEOS.py`).

Смесь инвариантных/known-answer проверок (детерминизм, физичность выбора корня,
предел идеального газа, корректность аналитической dZ/dp) и характеризационных
пинов текущего численного вывода (`pytest.approx(rel=1e-6)`) на реальном составе
`KRSNL_PVTSIM`. Как и `tests/regression/test_flash.py`, характеризационные числа
не сверены с внешним эталоном — их цель ловить непреднамеренные изменения при
рефакторинге. Если число законно меняется — пересчитать пин осознанно.

Паттерн использования EOS (как во flash-тесте): выставить `composition.T = t_K`
до конструирования (пересчитывает температурно-зависимые BRS-параметры a(T)),
затем `BrusilovskiyEOS(composition, p_bar, t_K)` с тем же t.
"""

import numpy as np
import pytest

from calc_core.EOS.BrusilovskiyEOS import BrusilovskiyEOS

REL = 1e-6

# Точка одиночной (газовой) фазы — та же (P, T), что в flash-тесте
# `test_krsnl_single_phase_above_p_sat`, только T в Кельвинах напрямую.
P_SINGLE_BAR = 250.0
T_SINGLE_K = 390 + 273.15


def _eos_at(composition, p_bar, t_k):
    """Готовит состав к температуре и строит EOS (общий паттерн для тестов)."""
    composition.T = t_k
    return BrusilovskiyEOS(composition, p_bar, t_k)


def test_calc_eos_deterministic(krsnl_composition):
    """Два вызова `calc_eos()` подряд дают идентичный Z и набор корней (нет скрытого состояния/RNG)."""
    eos = _eos_at(krsnl_composition, P_SINGLE_BAR, T_SINGLE_K)

    z1, _ = eos.calc_eos()
    roots1 = eos.real_roots_eos.copy()
    z2, _ = eos.calc_eos()
    roots2 = eos.real_roots_eos.copy()

    assert z1 == z2
    assert np.array_equal(roots1, roots2)


def test_chosen_root_is_real_and_positive(krsnl_composition):
    """Выбранный физический корень присутствует среди действительных корней кубики и положителен."""
    eos = _eos_at(krsnl_composition, P_SINGLE_BAR, T_SINGLE_K)
    eos.calc_eos()

    assert eos.choosen_eos_root in list(eos.real_roots_eos)
    assert eos.z > 0.0


def test_ideal_gas_limit(krsnl_composition):
    """При очень низком давлении Z-фактор стремится к 1 (кубика вырождается в Z**3 - Z**2 = 0)."""
    eos = _eos_at(krsnl_composition, 0.01, T_SINGLE_K)
    z, _ = eos.calc_eos()

    assert z == pytest.approx(1.0, abs=1e-3)


def test_dz_dp_matches_finite_difference(krsnl_composition):
    """Аналитическая dZ/dp (`_calc_dz_dp`) совпадает с независимой центральной разностью Z(P).

    Прямо страхует исторически багованный участок производных (см. докстринг
    `_calc_dlogphi_dp_vector`): `_calc_dz_dp` — единственная аналитическая
    производная EOS по давлению, и она должна совпадать с численной.
    """
    eos = _eos_at(krsnl_composition, P_SINGLE_BAR, T_SINGLE_K)
    z, _ = eos.calc_eos()
    dzdp_analytic = eos._calc_dz_dp()

    h = max(abs(P_SINGLE_BAR) * 1e-5, 1e-6)

    def z_nearest(p_shifted):
        e = _eos_at(krsnl_composition, p_shifted, T_SINGLE_K)
        e.calc_eos()
        roots = e.real_roots_eos
        return float(roots[np.argmin(np.abs(roots - z))])

    dzdp_fd = (z_nearest(P_SINGLE_BAR + h) - z_nearest(P_SINGLE_BAR - h)) / (2.0 * h)

    assert dzdp_analytic == pytest.approx(dzdp_fd, rel=1e-5)


def test_eos_characterization_krsnl(krsnl_composition):
    """Характеризация: пин Z-фактора и ln(phi) первых трёх компонент (N2/CO2/C1) при (250 бар, 663.15 K)."""
    eos = _eos_at(krsnl_composition, P_SINGLE_BAR, T_SINGLE_K)
    eos.calc_eos()

    assert eos.z == pytest.approx(1.0574203409162748, rel=REL)

    ln_phi = eos.get_fugacity_coef_vector_by_root(eos.choosen_eos_root)
    assert float(ln_phi[0]) == pytest.approx(0.5761355819270607, rel=REL)  # N2
    assert float(ln_phi[1]) == pytest.approx(0.2121850114709165, rel=REL)  # CO2
    assert float(ln_phi[2]) == pytest.approx(0.3467129511030398, rel=REL)  # C1


def test_constructor_does_not_mutate_composition(krsnl_composition):
    """Конструирование EOS и `calc_eos()` не меняют переданный `Composition` (T и мольные доли)."""
    krsnl_composition.T = T_SINGLE_K
    t_before = krsnl_composition.T
    composition_before = dict(krsnl_composition.composition)

    eos = BrusilovskiyEOS(krsnl_composition, P_SINGLE_BAR, T_SINGLE_K)
    eos.calc_eos()

    assert krsnl_composition.T == t_before
    assert dict(krsnl_composition.composition) == composition_before
