"""Тесты решателя двухфазного равновесия (`calc_core/VLE/PhaseEquilibriumNewton.py`).

Фокус — на сигнале несходимости: раньше при исчерпании лимита итераций внешний
цикл Ньютона по ln(K) тихо возвращал последний (несошедшийся) результат; теперь
бросает `ConvergenceError`. Плюс проверка штатной сходимости на известной
двухфазной точке.

Температура в K напрямую; `composition.T = t_k` до конструирования (как в других
unit-тестах и во flash-тесте).
"""

import pytest

from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest
from calc_core.VLE.PhaseEquilibriumNewton import PhaseEquilibriumNewton
from calc_core.Utils.Errors import ConvergenceError

# Известная двухфазная точка KRSNL (та же, что во flash-тесте): P=100 бар, T=110°C.
P_TWO_PHASE_BAR = 100.0
T_RES_K = 110 + 273.15
EXPECTED_KEYS = {"yi_v", "xi_l", "Ki", "Fv", "Fl", "Z_v", "Z_l"}


def _make_solver(composition):
    """Готовит K-значения тестом стабильности и строит решатель равновесия (как это делает Flash)."""
    composition.T = T_RES_K
    stability = TwoPhaseStabilityTest(composition, P_TWO_PHASE_BAR, T_RES_K)
    stability.calculate_phase_stability()
    assert stability.stable is False  # предпосылка теста: точка действительно двухфазная
    return PhaseEquilibriumNewton(composition, P_TWO_PHASE_BAR, T_RES_K, stability.k_values_for_flash)


def test_converges_on_two_phase_point(krsnl_composition):
    """Штатная сходимость: возвращается полный dict, Fv в (0, 1), флаг `convergence` выставлен."""
    solver = _make_solver(krsnl_composition)

    result = solver.find_solve_loop()

    assert set(result.keys()) == EXPECTED_KEYS
    assert 0.0 < result["Fv"] < 1.0
    assert solver.convergence is True


def test_raises_convergence_error_when_iterations_exhausted(krsnl_composition):
    """При исчерпании лимита итераций поднимается `ConvergenceError`, а не тихий возврат мусора."""
    solver = _make_solver(krsnl_composition)
    solver._FUG_MAX_ITER = 0  # форсируем несходимость на первой же проверке лимита

    with pytest.raises(ConvergenceError):
        solver.find_solve_loop()
