"""
Фактические константы движка для read-only панели Settings.

Собирает актуальные значения констант, стандартных условий и критериев
сходимости из движка (`calc_core.Utils.Constants`, `Conditions.StandardConditions`)
Панель намеренно не сохраняет «настройки»: до появления явного EngineConfig
любые редактируемые поля создавали бы ложное впечатление, что расчёт изменён.
"""

from calc_core.Utils.EngineConfig import EngineConfig

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
    {
        "group": "Iteration limits and parallelism",
        "fields": [
            {"key": "STABILITY_MAX_ITERATIONS", "label": "Stability max iterations", "fmt": "%d"},
            {"key": "FLASH_RR_MAX_ITERATIONS", "label": "Flash Rachford-Rice max iterations", "fmt": "%d"},
            {"key": "FLASH_FUGACITY_MAX_ITERATIONS", "label": "Flash fugacity max iterations", "fmt": "%d"},
            {"key": "PARALLEL_JOBS", "label": "Parallel jobs", "fmt": "%d"},
        ],
    },
]

# все ключи схемы по порядку (для сериализации/сравнения)
ALL_KEYS = [f["key"] for grp in SCHEMA for f in grp["fields"]]


def engine_defaults() -> dict:
    """Актуальные значения из движка (то, чем он реально считает сейчас)."""
    config = EngineConfig.defaults()
    return {
        "CONSTANT_R": config.constant_r,
        _STD_P: config.standard_pressure_bar,
        _STD_T: config.standard_temperature_k,
        "TOL_TWO_PHASE_STABILITY_CONVERGENCE": config.stability_convergence_tolerance,
        "TOL_TWO_PHASE_STABILITY_CONVERGENCE_TRIVIAL_SOLUTION": config.stability_trivial_tolerance,
        "TOL_TWO_PHASE_FLASH_BISECTION_CONVERGENCE": config.flash_bisection_tolerance,
        "TOL_TWO_PHASE_FLASH_CONVERGENCE": config.flash_convergence_tolerance,
        "TOL_TWO_PHASE_FLASH_TRIVIAL_SOLUTION": config.flash_trivial_tolerance,
        "TOL_SAT_PRESSURE": config.saturation_pressure_tolerance,
        "STABILITY_MAX_ITERATIONS": config.stability_max_iterations,
        "FLASH_RR_MAX_ITERATIONS": config.flash_rr_max_iterations,
        "FLASH_FUGACITY_MAX_ITERATIONS": config.flash_fugacity_max_iterations,
        "PARALLEL_JOBS": config.parallel_jobs,
    }
