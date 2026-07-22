"""Кооперативная отмена в основных расчётных точках."""

import pytest

from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
from calc_core.PhaseEnvelope.PhaseEnvelopeGrid import PhaseEnvelopeGrid
from calc_core.PhaseEnvelope.PhaseEnvelopeSuccessiveSubstitution import PhaseEnvelopeSSM
from calc_core.Utils.Cancellation import CalculationCancelled, CancellationToken


def test_cancelled_model_flash_stops_before_mutating_work(krsnl_composition):
    model = CompositionalModel(krsnl_composition)
    token = CancellationToken()
    token.cancel()

    with pytest.raises(CalculationCancelled):
        model.flash(50.0, 20.0, cancellation_token=token)


def test_cancelled_ssm_stops_at_temperature_boundary(krsnl_composition):
    calc = PhaseEnvelopeSSM(krsnl_composition, 60.0, 80.0, 20.0, 400.0)
    token = CancellationToken()
    token.cancel()

    with pytest.raises(CalculationCancelled):
        calc.calculate(token)


def test_cancelled_grid_stops_before_first_point(krsnl_composition):
    calc = PhaseEnvelopeGrid(
        krsnl_composition, max_pressure=10.0, max_temperature=20.0,
        pressure_points=2, temperature_points=2,
    )
    token = CancellationToken()
    token.cancel()

    with pytest.raises(CalculationCancelled):
        calc.run(token)
