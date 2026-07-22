"""Физические проверки уже полученного результата расчёта.

Диагностики не заменяют ``ConvergenceError``: если численный метод не
сошёлся, результата нет. Этот модуль отмечает физически подозрительные
значения в формально полученном результате, чтобы вызывающий код мог показать
состояние «получено, но требует проверки».
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Iterable, Mapping

if TYPE_CHECKING:
    from calc_core.VLE.FlashResult import FlashResult, PhaseState


@dataclass(frozen=True)
class ResultWarning:
    """Одно машиночитаемое предупреждение о качестве результата."""

    code: str
    message: str
    field: str | None = None


@dataclass(frozen=True)
class ResultDiagnostics:
    """Набор предупреждений; отсутствие предупреждений означает ``ok``."""

    warnings: tuple[ResultWarning, ...] = ()

    @property
    def status(self) -> str:
        return "warning" if self.warnings else "ok"

    def to_dict(self) -> dict[str, Any]:
        return {
            "status": self.status,
            "warnings": [
                {"code": item.code, "message": item.message, "field": item.field}
                for item in self.warnings
            ],
        }

    @classmethod
    def from_dict(cls, value: object) -> ResultDiagnostics:
        if not isinstance(value, dict):
            return cls()
        raw = value.get("warnings")
        if not isinstance(raw, list):
            return cls()
        warnings = []
        for item in raw:
            if not isinstance(item, dict):
                continue
            code, message = item.get("code"), item.get("message")
            field = item.get("field")
            if isinstance(code, str) and isinstance(message, str):
                warnings.append(ResultWarning(
                    code, message, field if isinstance(field, str) else None,
                ))
        return cls(tuple(warnings))


def _finite(value: object) -> bool:
    return (isinstance(value, (int, float)) and not isinstance(value, bool)
            and math.isfinite(float(value)))


def _warning(code: str, message: str, field: str | None = None) -> ResultWarning:
    return ResultWarning(code=code, message=message, field=field)


def _check_composition(
    composition: object,
    field: str,
    warnings: list[ResultWarning],
    *,
    tolerance: float,
) -> None:
    if not isinstance(composition, Mapping) or not composition:
        warnings.append(_warning(
            "composition_missing", f"Состав {field} отсутствует.", field,
        ))
        return
    values = list(composition.values())
    if any(not _finite(value) for value in values):
        warnings.append(_warning(
            "non_finite", f"Состав {field} содержит нечисловые/бесконечные значения.", field,
        ))
        return
    if any(float(value) < 0.0 or float(value) > 1.0 for value in values):
        warnings.append(_warning(
            "composition_bounds", f"Доли компонентов {field} выходят за [0, 1].", field,
        ))
    total = sum(float(value) for value in values)
    if not math.isclose(total, 1.0, abs_tol=tolerance, rel_tol=0.0):
        warnings.append(_warning(
            "composition_sum", f"Сумма долей {field} равна {total:.8g}, а не 1.", field,
        ))


def _check_phase_properties(
    phase: PhaseState,
    field: str,
    warnings: list[ResultWarning],
) -> None:
    if not _finite(phase.mole_fraction) or float(phase.mole_fraction) <= 0.0:
        return
    for name, value in phase.properties.items():
        path = f"{field}.properties.{name}"
        if not _finite(value):
            warnings.append(_warning(
                "non_finite", f"Свойство {path} не является конечным числом.", path,
            ))
    for name in ("molar_volume", "molar_density", "density"):
        value = phase.properties.get(name)
        if _finite(value) and float(value) <= 0.0:
            path = f"{field}.properties.{name}"
            warnings.append(_warning(
                "non_positive_property", f"Свойство {path} должно быть больше 0.", path,
            ))


def diagnose_flash_result(
    result: FlashResult,
    feed_composition: Mapping[str, float],
    *,
    tolerance: float = 1e-6,
    extra_warnings: Iterable[ResultWarning] = (),
) -> ResultDiagnostics:
    """Проверяет конечность, доли фаз, составы и материальный баланс Flash."""
    warnings = list(extra_warnings)

    for field, value in (("pressure", result.pressure),
                         ("temperature", result.temperature)):
        if not _finite(value) or float(value) <= 0.0:
            warnings.append(_warning(
                "invalid_condition", f"{field} должно быть конечным и больше 0.", field,
            ))

    fractions = (result.vapor.mole_fraction, result.liquid.mole_fraction)
    for field, value in zip(("vapor.mole_fraction", "liquid.mole_fraction"), fractions):
        if not _finite(value):
            warnings.append(_warning(
                "non_finite", f"{field} не является конечным числом.", field,
            ))
        elif float(value) < 0.0 or float(value) > 1.0:
            warnings.append(_warning(
                "phase_fraction_bounds", f"{field} выходит за [0, 1].", field,
            ))
    if all(_finite(value) for value in fractions):
        total = sum(float(value) for value in fractions)
        if not math.isclose(total, 1.0, abs_tol=tolerance, rel_tol=0.0):
            warnings.append(_warning(
                "phase_fraction_sum", f"Сумма долей фаз равна {total:.8g}, а не 1.",
                "phase_fractions",
            ))

    _check_composition(result.vapor.composition, "vapor.composition", warnings,
                       tolerance=tolerance)
    _check_composition(result.liquid.composition, "liquid.composition", warnings,
                       tolerance=tolerance)
    _check_phase_properties(result.vapor, "vapor", warnings)
    _check_phase_properties(result.liquid, "liquid", warnings)

    yi, xi = result.vapor.composition, result.liquid.composition
    if (all(_finite(value) for value in fractions)
            and isinstance(yi, Mapping) and isinstance(xi, Mapping)):
        fv, fl = float(fractions[0]), float(fractions[1])
        residuals = []
        for component, z_value in feed_composition.items():
            y_value, x_value = yi.get(component), xi.get(component)
            if not all(_finite(value) for value in (z_value, y_value, x_value)):
                continue
            residuals.append(abs(float(z_value) - fv * float(y_value) - fl * float(x_value)))
        max_residual = max(residuals, default=0.0)
        if max_residual > tolerance:
            warnings.append(_warning(
                "material_balance",
                f"Максимальная невязка материального баланса равна {max_residual:.3e}.",
                "material_balance",
            ))

    return ResultDiagnostics(tuple(warnings))

