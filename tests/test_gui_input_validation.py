from gui.services import input_validation_service as validation


def test_flash_validation_rejects_nonphysical_pressure_and_temperature():
    errors = validation.validate_flash_inputs(0, -273.15)

    assert errors == {
        "P": "Pressure must be greater than 0 bar",
        "T": "Temperature must be above absolute zero",
    }


def test_experiment_validation_normalizes_values_and_reports_invalid_grid():
    params, errors = validation.validate_experiment_inputs(
        "separator", "400, 300\n200", 90, 450, "60, 50, 20")

    assert errors == {}
    assert params == {
        "pressures": [400.0, 300.0, 200.0],
        "T_c": 90.0,
        "P_res": 450.0,
        "stage_temps_c": [60.0, 50.0, 20.0],
    }

    _params, errors = validation.validate_experiment_inputs(
        "dle", "200, 300", 90, 0)
    assert errors["pressures"] == "Pressure stages must be strictly descending"
    assert errors["P_res"] == "P res must be greater than 0 bar"


def test_envelope_validation_checks_active_method_only():
    assert validation.validate_envelope_inputs("grid", {
        "grid_t_max_c": 250,
        "grid_p_max_bar": 700,
        "grid_t_points": 20,
        "grid_p_points": 30,
    }) == {}

    errors = validation.validate_envelope_inputs("ssm", {
        "t_min_c": 100,
        "t_max_c": 90,
        "t_step_c": 0,
        "p_max_bar": -1,
    })
    assert errors == {
        "t_max_c": "T max must be greater than T min",
        "t_step_c": "T step must be greater than 0",
        "p_max_bar": "P max must be greater than 0 bar",
    }
