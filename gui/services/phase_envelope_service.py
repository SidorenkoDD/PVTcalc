"""
Расчёт фазовой огибающей P-T (+ пластовое Psat) для GUI.

Обёртка над `CompositionalModel.phase_envelope` (`PhaseEnvelopeFacade`), который
сам делает независимую копию состава и логирует вызовы. Здесь: сборка
`CompositionalModel` из нормализованной глубокой копии состава, вызов
выбранного метода и приведение результата к простому JSON-совместимому словарю
(NaN → None).

Два метода (ключ `params["method"]`):
- `"ssm"` — огибающая как непрерывная кривая (марш по температуре,
  `PhaseEnvelopeSSM`); три ветки (bubble/dew_upper/dew_lower) сшиваются в одну
  линию без разрывов (`_joined_envelope`, как `PhaseEnvelopeSSM.plot`);
- `"grid"` — скан термодинамической стабильности по сетке P×T
  (`PhaseEnvelopeGrid`): каждая точка помечается однофазной/двухфазной.

Критическая точка (и любые особые точки — крикондентерма/криконденбар)
намеренно НЕ считаются и НЕ рисуются: их расчёт в движке сейчас некорректен
(см. заметку автора 2026-07-17) — вернёмся к ним отдельно.

Составляющие результата считаются независимо и «мягко» (несходимость одной
не роняет остальные):

- **огибающая** (`ssm`) — таблица по сетке температур; NaN на части строк
  больше не трактуется вслепую: `diagnostics.rows` отличает физически
  отсутствующую ветвь (`not_present`) от несходимости теста стабильности;
- **пластовое Psat** при пластовой температуре — интерполяцией по уже
  посчитанной SSM-огибающей (верхняя ветка bubble/dew_upper) в точке T_res;
  отдельный солвер не запускается, маркер ложится точно на нарисованную
  кривую. Если T_res вне посчитанного диапазона T — `None` (рисуется только
  вертикаль пластовой T). Никакого `new_methodv2` — весь расчёт идёт через SSM.

Не импортирует DearPyGui; оркестрация потока/прогресса/отмены — во View.
"""

import logging
import math

from calc_core.Composition.Composition import Composition
from calc_core.CompositionalModel.CompositionalModel import CompositionalModel
from calc_core.Utils.Cancellation import CancellationToken, ProgressCallback

logger = logging.getLogger(__name__)

# Колонки SSM-таблицы (см. PhaseEnvelopeSSM) и их подписи на графике/в таблице.
_CURVES = [
    ("Bubble_bar", "bubble", "Bubble point"),
    ("Dew_upper_bar", "dew_upper", "Dew point"),
    ("Dew_lower_bar", "dew_lower", "Dew point (retrograde)"),
]


def default_envelope_params(composition: Composition) -> dict:
    """
    Разумные параметры огибающей для первого прогона (обоих методов).

    Температурный диапазон центрируется относительно пластовой температуры
    (`composition.T`, K), чтобы двухфазная область попадала в окно. Содержит
    поля обоих методов (`ssm`/`grid`) — активный выбирается ключом `method`.
    """
    t_res_c = round(composition.T - 273.15, 2)
    return {
        "method": "ssm",
        # SSM (огибающая как кривая по маршу температуры)
        "t_min_c": -100.0,
        "t_max_c": float(max(200.0, round(t_res_c + 100.0))),
        "t_step_c": 10.0,
        "p_max_bar": 700.0,
        # Grid (скан стабильности P×T; T идёт от 0 °C до grid_t_max_c)
        "grid_t_max_c": float(max(300.0, round(t_res_c + 200.0))),
        "grid_p_max_bar": 700.0,
        "grid_t_points": 30,
        "grid_p_points": 30,
        "T_res_c": t_res_c,
    }


def _normalized_copy(composition: Composition) -> Composition:
    """Нормализованная независимая копия состава (общий объект не «уводится»)."""
    return composition.new_composition(composition.composition, deep_copy=True)


def _f(v):
    """float или None (NaN/нечисло → None) для JSON."""
    try:
        f = float(v)
    except (TypeError, ValueError):
        return None
    return None if f != f else f  # NaN != NaN


def _interp_reservoir_psat(temps: list, primary: list, t_res_c: float):
    """
    Давление насыщения при `t_res_c` (bar) линейной интерполяцией по верхней
    ветке SSM-огибающей (`primary` — bubble или dew_upper на каждой T).

    Возвращает `None`, если валидных точек < 2 или `t_res_c` вне их диапазона
    (за пределы не экстраполируем — маркер тогда не рисуется).
    """
    pts = sorted((t, p) for t, p in zip(temps, primary)
                 if t is not None and p is not None)
    if len(pts) < 2 or t_res_c < pts[0][0] or t_res_c > pts[-1][0]:
        return None
    for (t0, p0), (t1, p1) in zip(pts, pts[1:]):
        if t0 <= t_res_c <= t1:
            if t1 == t0:
                return p0
            return p0 + (p1 - p0) * (t_res_c - t0) / (t1 - t0)
    return None


def _joined_envelope(temps: list, raw: dict) -> dict:
    """
    Одна непрерывная линия огибающей (как `PhaseEnvelopeSSM.plot`): верхняя
    ветка (bubble, иначе dew_upper) по возрастанию T + нижняя ретроградная
    ветка (dew_lower) по убыванию T — три колонки сшиваются в один контур без
    разрывов.
    """
    upper = sorted((t, (b if b is not None else d))
                   for t, b, d in zip(temps, raw["bubble"], raw["dew_upper"])
                   if t is not None and (b is not None or d is not None))
    lower = sorted(((t, dl) for t, dl in zip(temps, raw["dew_lower"])
                    if t is not None and dl is not None), reverse=True)
    env_t = [t for t, _ in upper] + [t for t, _ in lower]
    env_p = [p for _, p in upper] + [p for _, p in lower]
    return {"T": env_t, "P": env_p}


def _reservoir(params: dict, psat) -> dict | None:
    """Маркер пластовой точки (`T_res_c` из параметров + переданное Psat)."""
    t_res_c = params.get("T_res_c")
    if t_res_c is None:
        return None
    return {"T_c": float(t_res_c), "P_sat": psat}


def run_envelope(
    composition: Composition,
    params: dict,
    *,
    cancellation_token: CancellationToken | None = None,
    progress_callback: ProgressCallback | None = None,
) -> dict:
    """
    Считает фазовую огибающую выбранным методом (`params["method"]`).

    - `"ssm"` (по умолчанию) — огибающая как непрерывная кривая (марш по T,
      `PhaseEnvelopeSSM`); результат содержит `envelope` (сшитая линия) +
      `curves` (ветки по отдельности, для таблицы) + `reservoir` (Psat
      интерполяцией по огибающей).
    - `"grid"` — скан термодинамической стабильности по сетке P×T
      (`PhaseEnvelopeGrid`); результат содержит `grid` (точки stable/unstable).

    Возвращает JSON-совместимый dict (числа float, NaN → None). Критическая
    точка/особые точки не считаются (см. докстринг модуля).
    """
    method = params.get("method", "ssm")
    if method == "grid":
        return _run_grid(composition, params, cancellation_token, progress_callback)
    return _run_ssm(composition, params, cancellation_token, progress_callback)


def _run_ssm(composition, params, cancellation_token=None, progress_callback=None) -> dict:
    """Огибающая методом SSM (непрерывная кривая) — см. `run_envelope`."""
    t_min_c = float(params["t_min_c"])
    t_max_c = float(params["t_max_c"])
    t_step_c = float(params["t_step_c"])
    p_max_bar = float(params["p_max_bar"])
    if not all(math.isfinite(v) for v in (t_min_c, t_max_c, t_step_c, p_max_bar)):
        raise ValueError("Envelope parameters must be finite")
    if t_min_c <= -273.15 or t_max_c <= t_min_c or t_step_c <= 0 or p_max_bar <= 0:
        raise ValueError("Invalid SSM temperature/pressure range")
    if (t_max_c - t_min_c) / t_step_c > 5000:
        raise ValueError("SSM grid is too large (maximum 5000 temperature steps)")

    comp = _normalized_copy(composition)
    model = CompositionalModel(comp)
    logger.info("Огибающая (SSM): T=[%s..%s] шаг %s °C, P_max=%s бар",
                t_min_c, t_max_c, t_step_c, p_max_bar)

    df = model.phase_envelope.ssm(
        t_min_c, t_max_c, t_step_c, p_max_bar, parallel=True,
        cancellation_token=cancellation_token,
        progress_callback=progress_callback,
    )
    diagnostics = df.attrs.get("diagnostics", {
        "method": "ssm",
        "partial": False,
        "failed_temperatures": [],
        "fallback": None,
        "rows": [],
    })

    temps = [_f(t) for t in df["Temp_C"].tolist()]
    raw: dict = {}   # ключ ветки -> список значений колонки, выровненный с temps
    curves: dict = {}
    for col, key, label in _CURVES:
        col_vals = [_f(v) for v in df[col].tolist()] if col in df.columns else [None] * len(temps)
        raw[key] = col_vals
        ts, ps = [], []
        for t, p in zip(temps, col_vals):
            if t is not None and p is not None:
                ts.append(t)
                ps.append(p)
        curves[key] = {"T": ts, "P": ps, "label": label}

    # таблица (полная сетка, NaN → None)
    columns = ["Temp_C"] + [c for c, _, _ in _CURVES if c in df.columns]
    rows = [[_f(df[c].iloc[i]) for c in columns] for i in range(len(df))]

    reservoir = None
    if params.get("T_res_c") is not None:
        primary = [b if b is not None else d
                   for b, d in zip(raw["bubble"], raw["dew_upper"])]
        reservoir = _reservoir(
            params, _interp_reservoir_psat(temps, primary, float(params["T_res_c"])))

    return {
        "method": "ssm",
        "envelope": _joined_envelope(temps, raw),
        "curves": curves,
        "diagnostics": diagnostics,
        "reservoir": reservoir,
        "table": {"columns": columns, "rows": rows},
        "t_range": [t_min_c, t_max_c],
    }


def _run_grid(composition, params, cancellation_token=None, progress_callback=None) -> dict:
    """Скан стабильности P×T (`PhaseEnvelopeGrid`) — см. `run_envelope`."""
    t_max_c = float(params.get("grid_t_max_c", 300.0))
    p_max_bar = float(params.get("grid_p_max_bar", 700.0))
    t_points = int(params.get("grid_t_points", 30))
    p_points = int(params.get("grid_p_points", 30))
    if not math.isfinite(t_max_c) or not math.isfinite(p_max_bar):
        raise ValueError("Grid parameters must be finite")
    if t_max_c <= 0 or p_max_bar <= 0 or t_points < 2 or p_points < 2:
        raise ValueError("Invalid phase-envelope grid range")
    if t_points > 200 or p_points > 200 or t_points * p_points > 40000:
        raise ValueError("Phase-envelope grid is too large (maximum 200x200)")

    comp = _normalized_copy(composition)
    model = CompositionalModel(comp)
    logger.info("Огибающая (grid): T=[0..%s] °C, P=[1..%s] бар, сетка %d×%d",
                t_max_c, p_max_bar, t_points, p_points)

    calc = model.phase_envelope.grid(
        max_temperature=t_max_c, max_pressure=p_max_bar,
        temperature_points=t_points, pressure_points=p_points,
        cancellation_token=cancellation_token,
        progress_callback=progress_callback,
    )

    T = [_f(t) for t in calc.result_temperature_arr]   # уже °C
    P = [_f(p) for p in calc.result_pressure_arr]       # бар
    flags = [_f(f) for f in calc.result_stability_flag_arr]  # 1.0 = двухфазно

    unstable: dict[str, list[float]] = {"T": [], "P": []}
    stable: dict[str, list[float]] = {"T": [], "P": []}
    rows = []
    for t, p, f in zip(T, P, flags):
        if t is None or p is None:
            continue
        two_phase = (f == 1.0)
        (unstable if two_phase else stable)["T"].append(t)
        (unstable if two_phase else stable)["P"].append(p)
        rows.append([t, p, two_phase])

    return {
        "method": "grid",
        "grid": {"unstable": unstable, "stable": stable},
        # для grid Psat не интерполируем (нет кривой) — только вертикаль T_res
        "reservoir": _reservoir(params, None),
        "table": {"columns": ["Temp_C", "Pressure_bar", "Two_phase"], "rows": rows},
        "t_range": [0.0, t_max_c],
    }
