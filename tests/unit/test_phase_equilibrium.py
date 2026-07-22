"""Тесты решателя двухфазного равновесия (`calc_core/VLE/PhaseEquilibriumNewton.py`).

Фокус — на сигнале несходимости: раньше при исчерпании лимита итераций внешний
цикл Ньютона по ln(K) тихо возвращал последний (несошедшийся) результат; теперь
бросает `ConvergenceError`. Плюс проверка штатной сходимости на известной
двухфазной точке.

Температура в K напрямую; `composition.T = t_k` до конструирования (как в других
unit-тестах и во flash-тесте).
"""

import math
from types import SimpleNamespace

import pytest

from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest
from calc_core.Utils.Errors import ConvergenceError, InputValidationError
from calc_core.VLE.PhaseEquilibriumNewton import PhaseEquilibriumNewton

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


def _binary_solver(k_values=None):
    composition = SimpleNamespace(composition={"A": 0.5, "B": 0.5})
    return PhaseEquilibriumNewton(
        composition,
        p=1.0,
        t=300.0,
        k_values=k_values or {"A": 2.0, "B": 0.5},
    )


def test_converges_on_two_phase_point(krsnl_composition):
    """Штатная сходимость: возвращается полный dict, Fv в (0, 1), флаг `convergence` выставлен."""
    solver = _make_solver(krsnl_composition)

    result = solver.find_solve_loop()

    assert set(result.keys()) == EXPECTED_KEYS
    assert 0.0 < result["Fv"] < 1.0
    assert solver.convergence is True
    assert result["Fv"] + result["Fl"] == pytest.approx(1.0)
    assert sum(result["xi_l"].values()) == pytest.approx(1.0)
    assert sum(result["yi_v"].values()) == pytest.approx(1.0)
    for component, feed_fraction in krsnl_composition.composition.items():
        reconstructed = (
            result["Fl"] * result["xi_l"][component]
            + result["Fv"] * result["yi_v"][component]
        )
        assert reconstructed == pytest.approx(feed_fraction, abs=1e-10)


def test_raises_convergence_error_when_iterations_exhausted(krsnl_composition):
    """При исчерпании лимита итераций поднимается `ConvergenceError`, а не тихий возврат мусора."""
    solver = _make_solver(krsnl_composition)
    solver._FUG_MAX_ITER = 0  # форсируем несходимость на первой же проверке лимита

    with pytest.raises(ConvergenceError):
        solver.find_solve_loop()


@pytest.mark.parametrize(
    ("p", "t", "k_values", "message"),
    [
        (0.0, 300.0, {"A": 2.0, "B": 0.5}, "p"),
        (1.0, 0.0, {"A": 2.0, "B": 0.5}, "t"),
        (1.0, 300.0, None, "не должны быть None"),
        (1.0, 300.0, {"A": 2.0}, "не совпадают"),
        (1.0, 300.0, {"A": 2.0, "B": 0.0}, "больше 0"),
        (1.0, 300.0, {"A": 2.0, "B": math.nan}, "конечным"),
    ],
)
def test_rejects_invalid_pressure_temperature_and_k_values(p, t, k_values, message):
    composition = SimpleNamespace(composition={"A": 0.5, "B": 0.5})

    with pytest.raises(InputValidationError, match=message):
        PhaseEquilibriumNewton(composition, p=p, t=t, k_values=k_values)


def test_rejects_k_values_without_two_phase_rachford_rice_bracket():
    solver = _binary_solver({"A": 2.0, "B": 1.5})

    with pytest.raises(InputValidationError, match="Рэчфорда-Райза"):
        solver.find_solve_newton()


@pytest.mark.parametrize(
    ("method_name", "message"),
    [
        ("find_solve_newton", r"Rachford-Rice \(Newton\)"),
        ("find_solve_bisection_v4", r"Rachford-Rice \(bisection\)"),
    ],
)
def test_rachford_rice_raises_when_its_own_iteration_limit_is_exhausted(
    method_name, message,
):
    solver = _binary_solver()
    solver._RR_MAX_ITER = 0

    with pytest.raises(ConvergenceError, match=message):
        getattr(solver, method_name)()


def test_trivial_solution_flag_distinguishes_k_near_one():
    near_one = _binary_solver({"A": math.exp(1e-6), "B": math.exp(-1e-6)})
    regular = _binary_solver()

    assert near_one.check_trivial_solution() is True
    assert near_one.trivial_solution is True
    assert regular.check_trivial_solution() is False
    assert regular.trivial_solution is False
