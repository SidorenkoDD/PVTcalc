"""Перекрёстная верификация канонических алгоритмов фазовой огибающей.

Значения KRSNL ниже — внутренний regression baseline до получения чистого
экспорта PVTSim (R0.3). Физическая проверка не полагается только на числа:
каждая одноточечная граница подтверждается сменой флага устойчивости, а SSM
и Newton сравниваются как независимые алгоритмы на общей сетке температур.
"""

import numpy as np
import pytest

from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
from calc_core.PhaseEnvelope.BubblePointPressure import BubblePointCalculator
from calc_core.PhaseEnvelope.DewPressure import DewPointCalculator
from calc_core.PhaseEnvelope.Diagnostics import (
    NOT_PRESENT,
    NOT_SEARCHED,
    OK,
    STABILITY_NONCONVERGENCE,
    VERIFICATION_FAILED,
)
from calc_core.PhaseEnvelope.PhaseEnvelopeGrid import PhaseEnvelopeGrid
from calc_core.PhaseEnvelope.PhaseEnvelopeNewton import (
    PhaseEnvelopeNewton,
    SaturationPointNewton,
    verify_saturation_boundary,
)
from calc_core.PhaseEnvelope.PhaseEnvelopeSuccessiveSubstitution import (
    PhaseEnvelopeSSM,
    SaturationPointSSM,
)
from calc_core.Utils.Errors import ConvergenceError, InputValidationError


def _copy(composition):
    return composition.new_composition(composition.composition, deep_copy=True)


def test_bubble_point_matches_baseline_and_is_a_stability_boundary(krsnl_composition):
    t_k = 110.0 + 273.15
    calc = BubblePointCalculator(_copy(krsnl_composition), t_k)

    pressure = calc.calculate()

    assert calc.converged is True
    assert pressure == pytest.approx(185.0525818700687, rel=1e-9)
    assert abs(calc.result_F) < 1e-7
    assert verify_saturation_boundary(calc.composition, t_k, pressure, -1.0) == "bubble"
    model = CompositionalModel(_copy(krsnl_composition))
    assert model.phase_envelope.bubble_point(t_k) == pytest.approx(pressure, rel=1e-12)


def test_lower_dew_matches_baseline_and_is_a_stability_boundary(krsnl_composition):
    t_k = 400.0 + 273.15
    calc = DewPointCalculator(_copy(krsnl_composition), t_k, dew_point_type="lower")

    pressure = calc.calculate()

    assert calc.converged is True
    assert pressure == pytest.approx(3.170747795320686, rel=1e-8)
    assert abs(calc.result_F) < 1e-8
    assert verify_saturation_boundary(calc.composition, t_k, pressure, 1.0) == "dew"
    model = CompositionalModel(_copy(krsnl_composition))
    assert model.phase_envelope.dew_point(
        t_k, dew_point_type="lower",
    ) == pytest.approx(pressure, rel=1e-12)


def test_false_upper_dew_root_is_rejected_by_verifier_and_facade(krsnl_composition):
    """`converged=True` у сырого Newton-решателя ещё не доказывает границу фаз."""
    t_k = 400.0 + 273.15
    raw = DewPointCalculator(_copy(krsnl_composition), t_k, dew_point_type="upper")
    false_root = raw.calculate()

    assert raw.converged is True
    assert false_root == pytest.approx(579.0845147502829, rel=1e-8)

    verifier = SaturationPointNewton(_copy(krsnl_composition), t_k, p_max_bar=800.0)
    assert verifier._verify(false_root, nudge_sign=-1.0) is None
    assert verifier.last_verification_status == VERIFICATION_FAILED

    model = CompositionalModel(_copy(krsnl_composition))
    with pytest.raises(ConvergenceError, match="непроверенный корень"):
        model.phase_envelope.dew_point(t_k, dew_point_type="upper")


def test_ssm_does_not_publish_midpoint_after_stability_failure(
    krsnl_composition, monkeypatch,
):
    """Две несошедшиеся проверки внутри бисекции не превращаются в fake Psat."""
    finder = SaturationPointSSM(
        _copy(krsnl_composition), 110.0 + 273.15,
        p_max_bar=2.0, p_min_bar=1.0,
    )

    def fake_stability(p):
        if p == pytest.approx(2.0):
            return True
        if p == pytest.approx(1.0):
            return False
        finder._stability_failed = True
        return None

    monkeypatch.setattr(finder, "_is_stable", fake_stability)
    result = finder._bisect_stability_boundary(2.0, 1.0, (2.0, 1.0))

    assert result is None
    assert finder.find_bubble_point() is None
    assert finder.last_primary_status == STABILITY_NONCONVERGENCE


def test_ssm_and_newton_agree_on_common_bubble_grid(krsnl_composition):
    kwargs = dict(
        t_min_c=60.0,
        t_max_c=140.0,
        t_step_c=20.0,
        p_max_bar=400.0,
    )
    ssm = PhaseEnvelopeSSM(_copy(krsnl_composition), **kwargs)
    newton = PhaseEnvelopeNewton(
        _copy(krsnl_composition), **kwargs, bootstrap_lower_dew=False,
    )

    ssm_df = ssm.calculate()
    newton_df = newton.calculate()

    assert newton_df["Temp_C"].tolist() == ssm_df["Temp_C"].tolist()
    np.testing.assert_allclose(
        newton_df["Bubble_bar"].to_numpy(dtype=float),
        ssm_df["Bubble_bar"].to_numpy(dtype=float),
        rtol=1e-3,
    )
    assert all(row["primary_status"] == OK for row in ssm_df.attrs["diagnostics"]["rows"])
    assert all(
        row["secondary_status"] == NOT_PRESENT
        for row in ssm_df.attrs["diagnostics"]["rows"]
    )
    assert ssm_df.attrs["diagnostics"]["partial"] is False
    assert all(row["primary_status"] == OK for row in newton_df.attrs["diagnostics"]["rows"])
    assert all(
        row["secondary_status"] == NOT_SEARCHED
        for row in newton_df.attrs["diagnostics"]["rows"]
    )
    assert newton_df.attrs["diagnostics"]["fallback"] == "ssm"

    parallel_df = PhaseEnvelopeNewton(
        _copy(krsnl_composition), **kwargs, bootstrap_lower_dew=False,
    ).calculate_parallel(n_jobs=1)
    np.testing.assert_allclose(
        parallel_df["Bubble_bar"].to_numpy(dtype=float),
        newton_df["Bubble_bar"].to_numpy(dtype=float),
        rtol=1e-12,
    )
    assert len(parallel_df.attrs["diagnostics"]["rows"]) == len(parallel_df)


def test_grid_scan_has_expected_order_flags_and_no_composition_mutation(krsnl_composition):
    composition = _copy(krsnl_composition)
    original_temperature = composition.T
    grid = PhaseEnvelopeGrid(
        composition,
        max_pressure=400.0,
        max_temperature=200.0,
        pressure_points=5,
        temperature_points=5,
    )

    grid.run_parallel(n_jobs=1)

    assert composition.T == original_temperature
    np.testing.assert_allclose(
        np.asarray(grid.result_pressure_arr).reshape(5, 5)[0],
        [1.0, 100.75, 200.5, 300.25, 400.0],
    )
    np.testing.assert_allclose(
        np.asarray(grid.result_temperature_arr).reshape(5, 5)[:, 0],
        [0.0, 50.0, 100.0, 150.0, 200.0],
    )
    assert np.asarray(grid.result_stability_flag_arr).reshape(5, 5).tolist() == [
        [1.0, 0.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 0.0, 0.0],
        [1.0, 1.0, 1.0, 0.0, 0.0],
    ]


def test_phase_envelope_facade_rejects_nonpositive_temperature_step(krsnl_composition):
    model = CompositionalModel(_copy(krsnl_composition))

    with pytest.raises(InputValidationError, match="t_step_c"):
        model.phase_envelope.ssm(
            t_min_c=60.0,
            t_max_c=140.0,
            t_step_c=0.0,
            p_max_bar=400.0,
            parallel=False,
        )
