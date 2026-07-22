"""Тесты единого публичного фасада CompositionalModel.flash."""

import copy

import pytest

from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
from calc_core.Utils.Errors import InputValidationError
from calc_core.VLE.FlashResult import FlashResult


def test_model_flash_returns_typed_two_phase_result_without_mutating_composition(
    krsnl_composition,
):
    before_temperature = krsnl_composition.T
    before_data = copy.deepcopy(krsnl_composition.composition_data)
    model = CompositionalModel(krsnl_composition)

    result = model.flash(50.0, 20.0)

    assert isinstance(result, FlashResult)
    assert result.is_two_phase is True
    assert result.temperature == pytest.approx(293.15)
    assert result.pressure == pytest.approx(50.0)
    assert result.quality_status == "ok"
    assert krsnl_composition.T == before_temperature
    assert krsnl_composition.composition_data == before_data
    history = model.result_store_object.get_by_module("Flash")
    assert len(history) == 1
    assert history[0].params == {"p_bar": 50.0, "t_celsius": 20.0}
    assert history[0].data is result


def test_model_flash_returns_typed_single_phase_result(krsnl_composition):
    result = CompositionalModel(krsnl_composition).flash(200.0, 80.0)

    assert isinstance(result, FlashResult)
    assert result.is_two_phase is False
    assert result.phase_type == "ambiguous"
    assert result.liquid.mole_fraction == pytest.approx(1.0)
    assert result.vapor.mole_fraction == pytest.approx(0.0)


@pytest.mark.parametrize(
    ("p_bar", "t_celsius"),
    [(0.0, 20.0), (50.0, -273.15)],
)
def test_model_flash_validates_engineering_units(
    krsnl_composition, p_bar, t_celsius,
):
    model = CompositionalModel(krsnl_composition)

    with pytest.raises(InputValidationError):
        model.flash(p_bar, t_celsius)
