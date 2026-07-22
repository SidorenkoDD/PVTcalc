"""Контракт immutable EngineConfig и его передачи в основной Flash-путь."""

import math
from types import SimpleNamespace

import numpy as np
import pytest

from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest
from calc_core.Utils.Conditions import StandardConditions
from calc_core.Utils.Constants import CONSTANT_R, TOL_TWO_PHASE_FLASH_CONVERGENCE
from calc_core.Utils.EngineConfig import EngineConfig
from calc_core.Utils.Errors import InputValidationError
from calc_core.VLE.PhaseEquilibriumNewton import PhaseEquilibriumNewton
from gui.services.flash_service import restore_flash_result, snapshot_flash_result


def test_defaults_mirror_existing_engine_values():
    config = EngineConfig.defaults()
    standard = StandardConditions()

    assert config.constant_r == CONSTANT_R
    assert config.flash_convergence_tolerance == TOL_TWO_PHASE_FLASH_CONVERGENCE
    assert config.standard_pressure_bar == standard.p
    assert config.standard_temperature_k == standard.t


def test_config_is_immutable_and_rejects_invalid_limits():
    config = EngineConfig.defaults()

    with pytest.raises((AttributeError, TypeError)):
        config.parallel_jobs = 1

    with pytest.raises(InputValidationError):
        EngineConfig(flash_rr_max_iterations=0)


def test_config_changes_concrete_solver_tolerance():
    composition = SimpleNamespace(composition={"A": 0.5, "B": 0.5})
    config = EngineConfig(flash_convergence_tolerance=3.0)
    solver = PhaseEquilibriumNewton(
        composition, p=1.0, t=300.0,
        k_values={"A": 2.0, "B": 0.5}, config=config,
    )
    solver._ri_arr = np.array([2.0, 1.0])

    assert solver.config is config
    assert solver.check_convergence_ri() is True
    assert solver._FUG_MAX_ITER == config.flash_fugacity_max_iterations


def test_stability_uses_configured_tolerance(krsnl_composition):
    krsnl_composition.T = 383.15
    config = EngineConfig(stability_convergence_tolerance=1.0)
    stability = TwoPhaseStabilityTest(
        krsnl_composition, p=100.0, t=383.15, config=config,
    )

    assert stability._check_convergence(np.array([1.5])) is True
    assert stability.config is config


def test_model_and_result_snapshot_keep_the_same_config(krsnl_composition):
    config = EngineConfig(parallel_jobs=1, flash_rr_max_iterations=321)
    result = CompositionalModel(krsnl_composition, config=config).flash(50.0, 20.0)

    assert result.engine_config is config
    snap = snapshot_flash_result(result)
    assert snap["engine_config"]["flash_rr_max_iterations"] == 321

    restored = restore_flash_result(snap)
    assert restored.engine_config == config
    assert math.isclose(
        restored.engine_config.constant_r, config.constant_r,
    )
