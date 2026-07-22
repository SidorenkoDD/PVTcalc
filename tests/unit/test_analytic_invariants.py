"""Аналитические инварианты EOS/VLE на синтетических данных.

Эти тесты намеренно не создают ``Composition`` и не читают ``DB.json`` или
``models.json``. Минимальные объекты содержат только данные, которые нужны
проверяемому алгоритму. Благодаря этому падение локализуется в алгебре
EOS/Rachford–Rice/идентификации фазы, а не в корреляциях реальных компонентов.
"""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from calc_core.Composition.Composition import _normalize_composition
from calc_core.EOS.BrusilovskiyEOS import BrusilovskiyEOS
from calc_core.PhaseStability.PhaseIdentificator import PhaseIdentificator
from calc_core.VLE.PhaseEquilibriumNewton import PhaseEquilibriumNewton


def _pure_synthetic_composition() -> SimpleNamespace:
    """Минимальный однокомпонентный вход EOS без компонентной БД.

    Значения a/b/c/d выбраны конечными и различными, чтобы формула летучести
    не вырождалась в деление 0/0. При стремящемся к нулю давлении любые такие
    конечные параметры должны давать предел идеального газа.
    """
    component = "X"
    return SimpleNamespace(
        composition={component: 1.0},
        composition_data={
            "a": {component: 2.0},
            "b": {component: 0.03},
            "c": {component: 0.12},
            "d": {component: 0.01},
            "peneloux_correction": {component: 0.0},
            "bip": {component: {component: 0.0}},
        },
    )


def _binary_rr_solver() -> PhaseEquilibriumNewton:
    """Бинарная задача с точным решением Rachford–Rice Fv = 1/2."""
    feed = SimpleNamespace(composition={"A": 0.5, "B": 0.5})
    return PhaseEquilibriumNewton(
        feed,
        p=1.0,
        t=300.0,
        k_values={"A": 2.0, "B": 0.5},
    )


def test_synthetic_normalization_has_exact_ratios() -> None:
    """Нормировка сохраняет отношение исходных количеств 2:3."""
    normalized = _normalize_composition({"A": 2.0, "B": 3.0})

    assert normalized == pytest.approx({"A": 0.4, "B": 0.6})
    assert sum(normalized.values()) == pytest.approx(1.0)


def test_pure_component_reaches_ideal_gas_limit() -> None:
    """При P → 0 однокомпонентная EOS даёт Z → 1, phi → 1 и f → P."""
    pressure = 1e-8
    eos = BrusilovskiyEOS(_pure_synthetic_composition(), pressure, 350.0)

    z, log_fugacity = eos.calc_eos()
    log_phi = eos.get_fugacity_coef_vector_by_root(z)

    assert z == pytest.approx(1.0, abs=1e-10)
    assert float(log_phi[0]) == pytest.approx(0.0, abs=1e-10)
    assert float(log_fugacity[0]) == pytest.approx(math.log(pressure), abs=1e-10)


def test_pure_component_selects_ideal_gas_root_by_minimum_gibbs_energy() -> None:
    """Из нескольких действительных корней выбирается физичный корень Z≈1."""
    eos = BrusilovskiyEOS(_pure_synthetic_composition(), 1e-8, 350.0)

    z, _ = eos.calc_eos()

    assert len(eos.real_roots_eos) == 3
    assert z == pytest.approx(1.0, abs=1e-10)
    assert eos.normalized_gibbs_energy[eos._chosen_root_index] == pytest.approx(
        np.min(eos.normalized_gibbs_energy)
    )


def test_binary_rachford_rice_has_exact_half_vapor_solution() -> None:
    """Для z=(1/2,1/2), K=(2,1/2) корень RR аналитически равен 1/2."""
    solver = _binary_rr_solver()

    vapor_fraction = solver.find_solve_newton()

    assert vapor_fraction == pytest.approx(0.5, abs=1e-12)
    assert solver._rr_sum(vapor_fraction) == pytest.approx(0.0, abs=1e-12)


def test_binary_phase_compositions_close_material_balance() -> None:
    """Точные x/y нормированы и восстанавливают feed: z = L*x + V*y."""
    solver = _binary_rr_solver()
    solver.fv = solver.find_solve_newton()
    solver.L = 1.0 - solver.fv

    liquid, vapor = solver.define_xi_l_yi_v()

    assert liquid == pytest.approx({"A": 1.0 / 3.0, "B": 2.0 / 3.0})
    assert vapor == pytest.approx({"A": 2.0 / 3.0, "B": 1.0 / 3.0})
    assert sum(liquid.values()) == pytest.approx(1.0)
    assert sum(vapor.values()) == pytest.approx(1.0)
    for component, feed_fraction in solver.zi.items():
        reconstructed = solver.L * liquid[component] + solver.fv * vapor[component]
        assert reconstructed == pytest.approx(feed_fraction, abs=1e-12)


@pytest.mark.parametrize(
    ("s_v", "s_l", "molar_volume", "expected"),
    [
        (0.9, 1.0, 20.0, "vapor"),
        (1.1, 1.0, 1.0, "liquid"),
        (0.9, 1.0, 1.0, "ambiguous"),
    ],
)
def test_single_phase_identification_combines_two_independent_heuristics(
    s_v: float,
    s_l: float,
    molar_volume: float,
    expected: str,
) -> None:
    """Идентификатор возвращает фазу только при согласии stability/volume."""
    stability = SimpleNamespace(
        stable=True,
        S_v=s_v,
        S_l=s_l,
        S_v_rounded=round(s_v, 6),
        S_l_rounded=round(s_l, 6),
        liquid_eos=SimpleNamespace(
            _b=np.array([1.0]),
            composition_vector=np.array([1.0]),
        ),
    )
    flash = SimpleNamespace(
        phase_stability_object=stability,
        one_phase_stability_props={"molar_volume": molar_volume},
    )

    assert PhaseIdentificator(flash).identify_phase() == expected
