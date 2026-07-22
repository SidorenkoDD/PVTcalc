"""Immutable configuration of numerical engine defaults.

The defaults intentionally mirror the values that used to be read directly
from ``Constants`` or hard-coded in the main Flash solvers.  Keeping the
snapshot JSON-shaped makes a calculation reproducible without making the GUI
settings editable prematurely.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass
from typing import Any, Mapping

from calc_core.Utils.Conditions import StandardConditions
from calc_core.Utils.Constants import (
    CONSTANT_R,
    TOL_SAT_PRESSURE,
    TOL_TWO_PHASE_FLASH_BISECTION_CONVERGENCE,
    TOL_TWO_PHASE_FLASH_CONVERGENCE,
    TOL_TWO_PHASE_FLASH_TRIVIAL_SOLUTION,
    TOL_TWO_PHASE_STABILITY_CONVERGENCE,
    TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION,
)
from calc_core.Utils.Errors import InputValidationError


@dataclass(frozen=True, slots=True)
class EngineConfig:
    """All numerical defaults needed to describe one engine run."""

    constant_r: float = CONSTANT_R
    standard_pressure_bar: float = 1.01325
    standard_temperature_k: float = 293.15

    stability_convergence_tolerance: float = TOL_TWO_PHASE_STABILITY_CONVERGENCE
    stability_trivial_tolerance: float = TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION
    stability_max_iterations: int = 100000

    flash_bisection_tolerance: float = TOL_TWO_PHASE_FLASH_BISECTION_CONVERGENCE
    flash_convergence_tolerance: float = TOL_TWO_PHASE_FLASH_CONVERGENCE
    flash_trivial_tolerance: float = TOL_TWO_PHASE_FLASH_TRIVIAL_SOLUTION
    flash_rr_epsilon: float = 1e-12
    flash_rr_newton_tolerance: float = 1e-10
    flash_rr_max_iterations: int = 1000
    flash_fugacity_max_iterations: int = 1000

    saturation_pressure_tolerance: float = TOL_SAT_PRESSURE
    parallel_jobs: int = -1

    def __post_init__(self) -> None:
        positive_float_fields = (
            "constant_r", "standard_pressure_bar", "standard_temperature_k",
            "stability_convergence_tolerance", "stability_trivial_tolerance",
            "flash_bisection_tolerance", "flash_convergence_tolerance",
            "flash_trivial_tolerance", "flash_rr_epsilon",
            "flash_rr_newton_tolerance", "saturation_pressure_tolerance",
        )
        for name in positive_float_fields:
            value = getattr(self, name)
            if isinstance(value, bool) or not math.isfinite(float(value)) or float(value) <= 0.0:
                raise InputValidationError(f"EngineConfig.{name} должно быть конечным и больше 0.")

        for name in ("stability_max_iterations", "flash_rr_max_iterations",
                     "flash_fugacity_max_iterations"):
            value = getattr(self, name)
            if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
                raise InputValidationError(f"EngineConfig.{name} должно быть целым числом больше 0.")

        if (isinstance(self.parallel_jobs, bool) or not isinstance(self.parallel_jobs, int)
                or self.parallel_jobs == 0 or self.parallel_jobs < -1):
            raise InputValidationError(
                "EngineConfig.parallel_jobs должно быть целым: -1 или положительное число."
            )

    @classmethod
    def defaults(cls) -> "EngineConfig":
        """Return a fresh immutable config with the current engine defaults."""
        standard = StandardConditions()
        return cls(
            standard_pressure_bar=standard.p,
            standard_temperature_k=standard.t,
        )

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-shaped version suitable for a result snapshot."""
        return {"schema_version": 1, **asdict(self)}

    @classmethod
    def from_dict(cls, value: object) -> "EngineConfig":
        """Restore a config, falling back to defaults for legacy snapshots."""
        if not isinstance(value, Mapping):
            return cls.defaults()
        data = dict(value)
        data.pop("schema_version", None)
        defaults = cls.defaults()
        allowed = set(defaults.to_dict()) - {"schema_version"}
        data = {key: data[key] for key in allowed if key in data}
        try:
            return cls(**data)
        except (TypeError, ValueError, InputValidationError):
            return defaults
