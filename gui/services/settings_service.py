"""
Фактические константы движка для read-only панели Settings.

Собирает актуальные значения констант, стандартных условий и критериев
сходимости из движка (`calc_core.Utils.Constants`, `Conditions.StandardConditions`)
Панель намеренно не сохраняет «настройки»: до появления явного EngineConfig
любые редактируемые поля создавали бы ложное впечатление, что расчёт изменён.
"""

from calc_core.Utils import Constants
from calc_core.Utils.Conditions import StandardConditions

# Ключи стандартных условий (не константы модуля — читаются из StandardConditions).
_STD_P = "STD_P"
_STD_T = "STD_T"

# Схема панели: группы -> поля (key, подпись, формат отображения).
SCHEMA: list[dict] = [
    {
        "group": "Constants",
        "fields": [
            {"key": "CONSTANT_R", "label": "Universal gas constant R", "fmt": "%.6g"},
        ],
    },
    {
        "group": "Standard conditions",
        "fields": [
            {"key": _STD_P, "label": "Standard pressure, bar", "fmt": "%.6g"},
            {"key": _STD_T, "label": "Standard temperature, K", "fmt": "%.6g"},
        ],
    },
    {
        "group": "Convergence criteria (tolerances)",
        "fields": [
            {"key": "TOL_TWO_PHASE_STABILITY_CONVERGENCE",
             "label": "Stability test convergence", "fmt": "%.3g"},
            {"key": "TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION",
             "label": "Stability trivial solution", "fmt": "%.3g"},
            {"key": "TOL_TWO_PHASE_FLASH_BISECTION_CONVERGENCE",
             "label": "Flash bisection (Rachford-Rice)", "fmt": "%.3g"},
            {"key": "TOL_TWO_PHASE_FLASH_CONVERGENCE",
             "label": "Flash fugacity convergence", "fmt": "%.3g"},
            {"key": "TOL_TWO_PHASE_FLASH_TRIVIAL_SOLUTION",
             "label": "Flash trivial solution", "fmt": "%.3g"},
            {"key": "TOL_SAT_PRESSURE",
             "label": "Saturation pressure convergence", "fmt": "%.3g"},
        ],
    },
]

# все ключи схемы по порядку (для сериализации/сравнения)
ALL_KEYS = [f["key"] for grp in SCHEMA for f in grp["fields"]]


def engine_defaults() -> dict:
    """Актуальные значения из движка (то, чем он реально считает сейчас)."""
    std = StandardConditions()
    values = {name: getattr(Constants, name)
              for name in ALL_KEYS if hasattr(Constants, name)}
    values[_STD_P] = std.p
    values[_STD_T] = std.t
    return values
