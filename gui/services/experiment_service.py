"""
Запуск стандартных PVT-экспериментов (CCE / DLE / Separator) для GUI.

Обёртка над `CompositionalModel.experiments` (`ExperimentsFacade`), который
сам делает независимую копию состава и логирует вызовы. Здесь: сборка
`CompositionalModel` из переданного состава (нормализованная глубокая копия,
чтобы общий объект не «уводился»), диспетчеризация по типу эксперимента и
приведение результата (`pandas.DataFrame`) к простому JSON-совместимому
словарю таблицы (числовые колонки; dict-колонки составов отбрасываются).

Не импортирует DearPyGui; оркестрация потока/прогресса/отмены — во View.
"""

import logging

import pandas as pd

from calc_core.Composition.Composition import Composition
from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
from calc_core.Utils.Validation import (
    validate_positive_pressure,
    validate_temperature_celsius,
)

logger = logging.getLogger(__name__)

# Метаданные типов экспериментов: подпись, нужные поля ввода, наборы колонок
# для таблицы «Main» и для графиков (существующие в результате отбираются).
EXPERIMENT_TYPES: dict[str, dict] = {
    "cce": {
        "label": "CCE",
        "title": "Constant Composition Expansion",
        "needs_p_res": False,
        "needs_stage_temps": False,
        "plot_y": ["v/v_sat", "Compressibility"],
        "main_columns": ["pressure", "v/v_sat", "v/v_res", "Compressibility",
                         "vapor_mole_frac", "liquid_mole_frac", "liquid_density"],
        "charts": ["v/v_sat", "v/v_res", "Compressibility", "liquid_density",
                   "vapor_density", "liquid_z"],
    },
    "dle": {
        "label": "DLE",
        "title": "Differential Liberation",
        "needs_p_res": True,
        "needs_stage_temps": False,
        "plot_y": ["Bo", "Rs"],
        "main_columns": ["pressure", "Bo", "Rs", "vapor_mole_frac",
                         "liquid_mole_frac", "liquid_density", "vapor_density"],
        "charts": ["Bo", "Rs", "liquid_density", "vapor_density",
                   "liquid_z", "vapor_z"],
    },
    "separator": {
        "label": "Separator",
        "title": "Multi-stage separation",
        "needs_p_res": True,
        "needs_stage_temps": True,
        "plot_y": ["Bo", "Rs"],
        "main_columns": ["pressure", "Bo", "Rs", "vapor_mole_frac",
                         "liquid_mole_frac", "liquid_density"],
        "charts": ["Bo", "Rs", "liquid_density", "vapor_density"],
    },
}

_X_COLUMN = "pressure"
# колонки, которые не имеет смысла строить как кривую
_NON_PLOT = {"pressure", "temperature", "is_two_phase"}


def default_pressures(composition: Composition) -> list[float]:
    """Разумная сетка давлений по убыванию (бар) для первого прогона."""
    return [400.0, 350.0, 300.0, 250.0, 200.0, 150.0, 100.0, 50.0]


def _normalized_copy(composition: Composition) -> Composition:
    """Нормализованная независимая копия состава (см. докстринг модуля)."""
    return composition.new_composition(composition.composition, deep_copy=True)


def run_experiment(composition: Composition, kind: str, params: dict) -> dict:
    """
    Считает эксперимент и возвращает таблицу как словарь.

    Parameters
    ----------
    kind : str
        Один из `EXPERIMENT_TYPES` (`cce`/`dle`/`separator`).
    params : dict
        `pressures` (list, бар), `T_c` (°C), для dle/separator — `P_res` (бар),
        для separator — `stage_temps_c` (list, °C, той же длины, что pressures).

    Returns
    -------
    dict
        `{"columns": [...], "rows": [[...]], "x": "pressure", "plot_y": [...]}`
        — только числовые колонки, NaN → None (JSON-совместимо).
    """
    if kind not in EXPERIMENT_TYPES:
        raise ValueError(f"Неизвестный эксперимент: {kind}")

    comp = _normalized_copy(composition)
    model = CompositionalModel(comp)
    pressures = [float(p) for p in params["pressures"]]
    t_c = float(params["T_c"])
    if len(pressures) < 2:
        raise ValueError("At least two pressure stages are required")
    validate_positive_pressure(pressures, name="pressures")
    if any(a <= b for a, b in zip(pressures, pressures[1:])):
        raise ValueError("Pressure stages must be strictly descending")
    validate_temperature_celsius(t_c, name="T_c")
    logger.info("Эксперимент %s: %d ступеней, T=%s°C", kind, len(pressures), t_c)

    if kind == "cce":
        # facade.cce ждёт температуру в K (см. его докстринг)
        df = model.experiments.cce(pressures, t_c + 273.15)
    elif kind == "dle":
        validate_positive_pressure(float(params["P_res"]), name="P_res")
        df = model.experiments.dle(pressures, float(params["P_res"]), t_c)
    else:  # separator
        stage_temps = [float(t) for t in params["stage_temps_c"]]
        if len(stage_temps) != len(pressures):
            raise ValueError("Stage temperature count must match pressure count")
        for stage_t in stage_temps:
            validate_temperature_celsius(stage_t, name="stage temperature")
        validate_positive_pressure(float(params["P_res"]), name="P_res")
        df = model.experiments.separator(pressures, stage_temps,
                                         float(params["P_res"]), t_c)

    return _dataframe_to_result(df, kind)


def _f(v):
    """float или None (NaN/нечисло → None) для JSON."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f  # NaN


def _comp_dict(d):
    """Словарь состава фазы {компонент: доля} → JSON-совместимый (или None)."""
    if not isinstance(d, dict):
        return None
    return {str(k): _f(v) for k, v in d.items()}


def _dataframe_to_result(df: pd.DataFrame, kind: str) -> dict:
    """
    DataFrame эксперимента → богатый словарь результата:
    числовая таблица (все скалярные колонки), наборы для «Main»/графиков,
    и посоставный разбор по ступеням (`stages`).
    """
    meta = EXPERIMENT_TYPES[kind]

    # посоставный разбор по ступеням (до отбрасывания dict-колонок)
    stages: list[dict] = []
    for _, r in df.iterrows():
        stages.append({
            "pressure": _f(r.get("pressure")),
            "temperature": _f(r.get("temperature")),
            "vapor": _comp_dict(r.get("vapor_composition")),
            "liquid": _comp_dict(r.get("liquid_composition")),
        })

    num_df = df.select_dtypes(include=["number", "bool"])
    ordered = ([_X_COLUMN] if _X_COLUMN in num_df.columns else [])
    ordered += [c for c in meta["plot_y"] if c in num_df.columns]
    ordered += [c for c in num_df.columns if c not in ordered]

    rows: list[list] = []
    for _, r in num_df[ordered].iterrows():
        row: list = []
        for c in ordered:
            v = r[c]
            if pd.isna(v):
                row.append(None)
            elif isinstance(v, bool) or str(type(v)).find("bool") >= 0:
                row.append(bool(v))
            else:
                row.append(float(v))
        rows.append(row)

    return {
        "kind": kind,
        "x": _X_COLUMN,
        "columns": ordered,
        "rows": rows,
        "main_columns": [c for c in meta["main_columns"] if c in ordered],
        "charts": [c for c in meta["charts"] if c in ordered],
        "plot_all": [c for c in ordered if c not in _NON_PLOT],
        "plot_y": [c for c in meta["plot_y"] if c in ordered],
        "stages": stages,
    }


def series_for_plot(result: dict, y_column: str) -> tuple[list, list]:
    """Пары (x, y) для графика по колонке `y_column`, отсортированные по x, без None."""
    cols = result["columns"]
    if result["x"] not in cols or y_column not in cols:
        return [], []
    xi = cols.index(result["x"])
    yi = cols.index(y_column)
    pairs = [(row[xi], row[yi]) for row in result["rows"]
             if row[xi] is not None and row[yi] is not None]
    pairs.sort(key=lambda p: p[0])
    return [p[0] for p in pairs], [p[1] for p in pairs]
