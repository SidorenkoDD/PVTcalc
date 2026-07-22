"""Контрольные значения и негативные сценарии корреляций C7+.

Все методы, зарегистрированные через ``get_correlation``/``get_required_params``,
проверяются на одном инженерно правдоподобном псевдокомпоненте:
M=200 g/mol, gamma=0.85, Tb=600 K. Это внутренний regression baseline формул
и единиц; независимую предметную сверку диапазонов применимости нужно проводить
по первичным публикациям перед добавлением новых методов в GUI.
"""

import math

import pytest

from calc_core.PlusComponents.AcentricFactor import AcentricFactorCorrelation
from calc_core.PlusComponents.CriticalPressure import CriticalPressureCorrelation
from calc_core.PlusComponents.CriticalTemperature import CriticalTemperatureCorrelation
from calc_core.PlusComponents.CriticalVolume import CriticalVolumeCorrelation
from calc_core.PlusComponents.kWatson import KWatsonCorrelation
from calc_core.PlusComponents.PlusComponentCorrelations import PlusComponentProperties
from calc_core.PlusComponents.ShiftParameter import ShiftParameterCorrelation
from calc_core.Utils.Errors import InputValidationError

BASE = {
    "gamma": 0.85,
    "Tb_K": 600.0,
    "M": 200.0,
    "Kw": 11.5,
    "Tc_K": 750.0,
    "Pc_bar": 20.0,
    "af": 0.6,
}


CASES = [
    # Critical temperature, K
    (CriticalTemperatureCorrelation, "roess", 774.3193121030324),
    (CriticalTemperatureCorrelation, "nokey", 776.0783965425120),
    (CriticalTemperatureCorrelation, "cavett", 777.8725925440609),
    (CriticalTemperatureCorrelation, "kesler lee", 768.6953271604939),
    (CriticalTemperatureCorrelation, "pedersen", 672.1458060261916),
    (CriticalTemperatureCorrelation, "standing", 692.1607730605405),
    (CriticalTemperatureCorrelation, "sim daubert", 770.4404187642627),
    (CriticalTemperatureCorrelation, "riazi daubert", 775.6777960979939),
    (CriticalTemperatureCorrelation, "mogoulas tassios", 698.1029411764706),
    (CriticalTemperatureCorrelation, "twu", 755.0121805692189),
    (CriticalTemperatureCorrelation, "watansiri owens starling", 786.6403197928239),
    # Critical pressure, bar
    (CriticalPressureCorrelation, "kesler lee", 14.485123999514302),
    (CriticalPressureCorrelation, "riazi daubert", 14.273061260280892),
    (CriticalPressureCorrelation, "cavett", 14.655568219497434),
    (CriticalPressureCorrelation, "pedersen", 18.824420694275585),
    (CriticalPressureCorrelation, "standing", 19.911751204415037),
    (CriticalPressureCorrelation, "sim daubert", 14.942350513525872),
    (CriticalPressureCorrelation, "mogoulas tassios", 19.709447488853726),
    # Acentric factor, dimensionless
    (AcentricFactorCorrelation, "kesler lee", 0.9579475000000008),
    (AcentricFactorCorrelation, "riazi al sahhaf", 0.6406689062551210),
    (AcentricFactorCorrelation, "edmister", 1.2203318257819404),
    (AcentricFactorCorrelation, "mogoulas tassios", 0.5324654853672584),
    # Critical volume, cm3/mol. Reid includes J -> bar*cm3 conversion.
    (CriticalVolumeCorrelation, "hall yarborough", 786.1578033624303),
    (CriticalVolumeCorrelation, "riazi daubert", 977.9010476741028),
    (CriticalVolumeCorrelation, "reid", 736.2040925243787),
    (CriticalVolumeCorrelation, "lohrenz", 817.4535580259999),
    # Watson factor and volume-shift, dimensionless
    (KWatsonCorrelation, "k watson", 12.070418447129624),
    (KWatsonCorrelation, "riazi daubert", 11.687124143753946),
    (ShiftParameterCorrelation, "jhaveri youngren", 0.12311920730819448),
    (ShiftParameterCorrelation, "srk", 0.2658579250000001),
    (ShiftParameterCorrelation, "pr", 0.14006674439999997),
    (ShiftParameterCorrelation, "", 0.0),
]


@pytest.mark.parametrize(
    ("correlation_class", "method", "expected"),
    CASES,
    ids=[f"{cls.__name__}-{method or 'zero'}" for cls, method, _ in CASES],
)
def test_every_registered_correlation_has_control_value(
    correlation_class, method, expected,
):
    correlation = correlation_class()
    required = correlation.get_required_params(method)
    kwargs = {name: BASE[name] for name in required}

    result = correlation.get_correlation(method.upper())(**kwargs)

    assert math.isfinite(result)
    assert result == pytest.approx(expected, rel=1e-12)


DEFAULT_CONFIG = {
    "critical_temperature": "pedersen",
    "critical_pressure": "riazi daubert",
    "acentric_factor": "riazi al sahhaf",
    "critical_volume": "hall yarborough",
    "Kw": "k watson",
    "shift_parameter": "jhaveri youngren",
}


def test_default_c7_pipeline_has_expected_units_and_values():
    properties = PlusComponentProperties(
        M=200.0, gamma=0.85, Tb=600.0,
        correlations_config=DEFAULT_CONFIG,
    )

    properties.calculate_all()

    assert properties.data == pytest.approx({
        "molar_mass": 200.0,
        "gamma": 0.85,
        "Tb": 600.0,
        "Kw": 12.070418447129624,
        "critical_temperature": 672.1458060261916,
        "critical_pressure": 14.273061260280892,
        "acentric_factor": 0.6406689062551210,
        "critical_volume": 786.1578033624303,
        "shift_parameter": 0.12311920730819448,
    })


def test_precomputed_properties_are_not_overwritten_by_calculate_all():
    properties = PlusComponentProperties(
        M=200.0,
        gamma=0.85,
        Tb=600.0,
        Tc=750.0,
        Pc=20.0,
        Kw=11.5,
        af=0.6,
        correlations_config={
            "critical_volume": "reid",
            "shift_parameter": "pr",
        },
    )

    properties.calculate_all()

    assert properties.data["critical_temperature"] == 750.0
    assert properties.data["critical_pressure"] == 20.0
    assert properties.data["Kw"] == 11.5
    assert properties.data["acentric_factor"] == 0.6
    assert properties.data["critical_volume"] == pytest.approx(736.2040925243787)
    assert properties.data["shift_parameter"] == pytest.approx(0.1400667444)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"M": 0.0, "gamma": 0.85, "Tb": 600.0}, "M"),
        ({"M": 200.0, "gamma": math.nan, "Tb": 600.0}, "gamma"),
        ({"M": 200.0, "gamma": 0.85, "Tb": -1.0}, "Tb"),
    ],
)
def test_aggregator_rejects_nonphysical_base_inputs(kwargs, message):
    with pytest.raises(InputValidationError, match=message):
        PlusComponentProperties(**kwargs, correlations_config=DEFAULT_CONFIG)


def test_aggregator_reports_unknown_property_method_and_missing_dependency():
    properties = PlusComponentProperties(
        M=200.0,
        gamma=0.85,
        Tb=600.0,
        correlations_config={"acentric_factor": "kesler lee"},
    )

    with pytest.raises(ValueError, match="Unknown property"):
        properties.calculate_property("not_a_property")
    with pytest.raises(ValueError, match="No correlation specified"):
        properties.calculate_property("critical_pressure")
    with pytest.raises(ValueError, match="Missing required parameter 'Pc_bar'"):
        properties.calculate_property("acentric_factor")


@pytest.mark.parametrize(
    "correlation",
    [
        CriticalTemperatureCorrelation(),
        CriticalPressureCorrelation(),
        AcentricFactorCorrelation(),
        CriticalVolumeCorrelation(),
        KWatsonCorrelation(),
        ShiftParameterCorrelation(),
    ],
    ids=lambda correlation: type(correlation).__name__,
)
def test_unknown_correlation_is_rejected(correlation):
    with pytest.raises(ValueError, match="Unknown correlation method"):
        correlation.get_correlation("not registered")
    with pytest.raises(ValueError, match="Unknown correlation method"):
        correlation.get_required_params("not registered")


def test_unimplemented_todo_correlations_are_not_publicly_registered():
    with pytest.raises(ValueError, match="Unknown correlation method"):
        CriticalPressureCorrelation().get_correlation("twu")
    with pytest.raises(ValueError, match="Unknown correlation method"):
        CriticalVolumeCorrelation().get_correlation("riedel")
    with pytest.raises(ValueError, match="Unknown correlation method"):
        ShiftParameterCorrelation().get_correlation("brs")


def test_standing_limits_raise_domain_specific_error_before_math_domain_error():
    with pytest.raises(ValueError, match="61.1"):
        CriticalPressureCorrelation().get_correlation("standing")(gamma=0.85, M=60.0)
    with pytest.raises(ValueError, match="71.2"):
        CriticalTemperatureCorrelation().get_correlation("standing")(gamma=0.85, M=71.2)
