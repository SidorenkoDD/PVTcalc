"""
Фазовая огибающая методом последовательных приближений (successive substitution).

Независимый порт алгоритма `PhaseDiagram_v4.py` на текущую архитектуру
(`calc_core.*`, `TwoPhaseStabilityTest` без параметра `eos`, температура через
`Composition.T`). Метод Ньютона по производным летучести (`BubblePointPressure.py`/
`DewPressure.py`, "Table 6.1/6.2" по мануалу PVTSim) хорошо сходится в отдельной
точке, но теряет устойчивость вблизи критической точки — там PVTSim рекомендует
именно метод последовательных приближений, реализованный здесь.

Как и в `PhaseDiagram_v4.py`, НАПРАВЛЕНИЕ поиска на каждой температуре
задаётся не физикой, а механикой бисекции: `find_bubble_point()` ищет от
`p_max` вниз, `find_dew_point()` — от `p_min` вверх, в границах уже найденной
верхней точки на этой же температуре (второй точки может не быть вовсе).
Однако, в отличие от первой версии этого модуля, физический ТИП найденной
точки (Bubble/Dew) определяется ОТДЕЛЬНО от направления поиска — критерием
Михельсена (сравнение `S_v`/`S_l`, см. `SaturationPointSSM.
_classify_point_type`), а не жёстко привязан к тому, каким поиском она
найдена. Это принципиально: за критической точкой состава (которую саму по
себе мы не считаем — см. CLAUDE.md) поиск "сверху вниз" физически находит уже
не bubble point, а верхнюю dew point — обе найденные на такой T точки
оказываются точками росы (верхняя и нижняя, ограничивающие ретроградную
двухфазную область), а не "bubble + dew". Обе ветки ищутся независимо на
каждом шаге именно поэтому — не потому что механически на T может быть ровно
одна bubble и одна dew, а потому что до критической точки состава существует
только верхняя точка (физически bubble), а после — до двух точек росы.

Отступления от буквальной логики `PhaseDiagram_v4.py` (оба — ради устойчивости,
которая явно требовалась в задаче):

1. Направление бисекции определяется ДИСКРЕТНЫМ флагом стабильности
   (`TwoPhaseStabilityTest.stable`), а не знаком `(Σy - 1)` из SS-уточнения.
   Изначальная реализация на знаке `(Σy-1)` работала нестабильно именно там,
   где она нужнее всего — вблизи критической/ретроградной области непрерывный
   SS-сигнал может быть немонотонным (несколько корней уравнения Σy=1 рядом
   друг с другом), и бисекция по его знаку либо сходится не к тому корню, либо
   "проскакивает" мимо настоящей второй (нижней) точки насыщения. Дискретный
   флаг "стабильно/нестабильно" такой неоднозначности не подвержен: пока на
   границах интервала флаги гарантированно разные, бисекция по флагу сходится
   к единственному разделяющему их переходу.
2. Continuation по температуре (давление, найденное на предыдущем шаге,
   используется как узкая стартовая попытка bracket'а для следующего — с
   безопасным откатом на полный `[p_min_bar, p_max_bar]`, если этот узкий
   bracket не содержит перехода стабильности). В `PhaseDiagram_v4.py` такой
   передачи между температурами не было вовсе (каждый шаг стартовал заново).

Оптимизации скорости (все — только ускоряют поиск, ни одна не ослабляет
гарантию корректности: полный safety-net — скан + бисекция по флагу — всегда
остаётся откатом, если ускоряющая попытка не подошла):

- `_wilson_bubble_pressure_estimate` — дешёвая (без единого вызова
  `TwoPhaseStabilityTest`) оценка по корреляции Вильсона используется как
  guess по умолчанию, если явный не передан (в т.ч. в параллельном расчёте,
  где нет continuation) — может быть грубой (проверено: для тяжёлых составов
  расхождение с истинным Bubble доходит до 5-7x), но это не мешает
  корректности — `quick_bracket` сам её отбраковывает при неудаче.
- `_find_transition_scanning` использует геометрически РАСТУЩИЙ шаг по
  log(P) (не равномерную сетку) — мелкое разрешение у стартовой точки
  (`p_from`), кратно (не построчно) растущее дальше. Дёшево именно в
  наиболее частом случае — когда второй точки насыщения на этой T нет вовсе
  (типично для большей части диапазона T) — вместо полного прохода по всей
  равномерной сетке распознаёт "перехода нет" за O(log(диапазон)) проверок.
- Допуск бисекции — ОТНОСИТЕЛЬНЫЙ (`bracket_tol_rel`, см. `SaturationPointSSM.
  __init__`), а не абсолютный. Обнаружено профилированием: вблизи самой
  границы фазового перехода `TwoPhaseStabilityTest` сходится экспоненциально
  медленнее ("критическое замедление"), и последние несколько итераций
  бисекции при слишком жёстком абсолютном допуске могут стоить ДЕСЯТКИ СЕКУНД
  на одну точку ради точности, которая на порядок превосходит практическую
  потребность (реальное расхождение с PVTSim на этом составе — ~0.03%, см.
  ноутбук). Это дало основной прирост скорости во всей серии оптимизаций.
- `calculate_parallel` делит диапазон T на непрерывные КУСКИ (по числу
  доступных воркеров), а не на отдельную задачу под каждую точку —
  continuation работает внутри куска (как в `calculate()`), "холодный" поиск
  без guess'а нужен только для первой точки каждого куска.
"""

import logging
from typing import NamedTuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, effective_n_jobs
from calc_core.Utils.Cancellation import CancellationToken, ProgressCallback, report_progress

from calc_core.Composition.Composition import Composition
from calc_core.PhaseEnvelope.Diagnostics import (
    NOT_PRESENT,
    OK,
    STABILITY_NONCONVERGENCE,
    EnvelopeDiagnostic,
    attach_diagnostics,
)
from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest
from calc_core.Utils.Errors import StopIterationError

logger = logging.getLogger(__name__)


class SaturationPoint(NamedTuple):
    """
    Найденная точка насыщения + её ФИЗИЧЕСКИЙ тип (см. `SaturationPointSSM.
    _classify_point_type`) — не путать с направлением поиска, которым точка
    была найдена (`find_bubble_point`/`find_dew_point`): за критической
    точкой состава направление "сверху вниз" находит уже верхнюю Dew, а не
    Bubble (см. докстринг модуля).
    """
    pressure: float
    point_type: str  # 'bubble' | 'dew'


def _wilson_bubble_pressure_estimate(composition: Composition, t_k: float) -> float:
    """
    Дешёвая аналитическая оценка давления насыщения по корреляции Вильсона —
    тот же приём, что и в `BubblePointPressure.py`/`DewPressure.py` для
    начального приближения перед Ньютоном: `P = Σ z_i Pc_i exp(5.373(1+ω_i)(1-Tc_i/T))`.

    Ни одного вызова `TwoPhaseStabilityTest` не требует — считается напрямую
    по `composition_data`. Используется только чтобы сориентировать
    бисекцию (`quick_bracket`), не как решение: оценка может быть заметно не
    точной вблизи критической точки/для газоконденсатных составов, но за
    корректность в этом случае отвечает scan-fallback в
    `_bisect_stability_boundary` — потеря точности этой оценки не может
    привести к неверному результату, только к более медленному его поиску.
    """
    z = composition.composition
    comp_data = composition.composition_data
    components = tuple(z.keys())

    z_arr = np.array([z[c] for c in components])
    pc = np.array([comp_data['critical_pressure'][c] for c in components])
    tc = np.array([comp_data['critical_temperature'][c] for c in components])
    omega = np.array([comp_data['acentric_factor'][c] for c in components])

    exp_term = np.exp(5.373 * (1.0 + omega) * (1.0 - tc / t_k))
    return float(np.sum(z_arr * pc * exp_term))


class SaturationPointSSM:
    """
    Поиск давления насыщения (Bubble) и давления начала конденсации (Dew)
    при фиксированной температуре методом последовательных приближений —
    бисекцией по дискретному флагу теста стабильности Михельсена.

    Порт `PhaseDiagram_v4.py::SaturationPressure` на `TwoPhaseStabilityTest
    (composition, p, t)` без параметра `eos`.
    """

    def __init__(
        self,
        composition: Composition,
        t_k: float,
        p_max_bar: float,
        p_min_bar: float = 0.1,
        p_guess_bubble: float | None = None,
        p_guess_dew: float | None = None,
        bracket_tol_rel: float = 1e-3,
        bracket_tol_abs: float = 1e-3,
        max_iter: int = 100,
    ):
        self.composition = composition
        self.composition.T = t_k
        self.t_k = t_k

        self.p_min_bar = p_min_bar
        self.p_max_bar = p_max_bar
        # Точность бисекции — ОТНОСИТЕЛЬНАЯ (`bracket_tol_rel`, доля от P), не
        # абсолютная: диапазон давлений здесь может охватывать 3-4 порядка
        # (от долей бара для нижней ретроградной точки до тысяч бар для
        # верхней) — единый абсолютный допуск либо избыточно точен на
        # верхнем конце (тратит лишние итерации бисекции без практической
        # пользы: перекрёстная сверка с PVTSim на этом составе и так сходится
        # только с точностью ~0.05 бар, см. ноутбук), либо слишком грубый на
        # нижнем. `bracket_tol_abs` — нижний пол на случай очень малых P
        # (вблизи `p_min_bar`), чтобы не требовать бессмысленно малую
        # абсолютную точность там, где относительная теряет смысл.
        #
        # Дефолт 1e-3 (не 1e-4) — НЕ произвольный выбор ради скорости: вблизи
        # самой границы фазового перехода тест стабильности Михельсена
        # сходится экспоненциально медленнее ("критическое замедление" — S_v/
        # S_l приближаются к 1.0 очень медленно чем ближе P к истинному
        # давлению насыщения). На этом составе замерено: последние несколько
        # итераций бисекции при tol=1e-4 добавляли до 50 СЕКУНД на одну точку
        # (T=400°C, самая тяжёлая точка диапазона — вблизи пика Bubble-кривой,
        # т.е. ближе всего к критической области). tol=1e-3 даёт почти весь
        # выигрыш в скорости (~2x на последовательном марше, ~3x на
        # параллельном относительно tol=1e-4) и при этом СОГЛАСОВАН между
        # последовательным и параллельным расчётом (разными путями бисекции)
        # с точностью 1e-3 отн. — проверено на полном диапазоне T на двух
        # составах. Более агрессивный tol=1e-2 даёт ещё на порядок скорости
        # (~7x на последовательном, ~10x на параллельном относительно 1e-4),
        # но вблизи пика/приближения к критической области (T≈220-260°C на
        # этом составе) там, где S(P) физически "плоская", результат заметно
        # зависит от конкретного пути бисекции (до ~0.6-0.7% расхождения
        # между sequential и parallel) — не баг, а отражение реальной
        # near-critical неопределённости, но это уже разброс того же порядка,
        # что и расхождение с PVTSim (~0.03%), поэтому не взят как дефолт;
        # можно передать `bracket_tol_rel=1e-2` явно, если вне критической
        # области и нужна максимальная скорость.
        self.bracket_tol_rel = bracket_tol_rel
        self.bracket_tol_abs = bracket_tol_abs
        self.max_iter = max_iter

        # Если явный guess (continuation с предыдущего шага по T) не задан —
        # используем дешёвую оценку Вильсона вместо полного скана "вслепую".
        # Она не обязана быть точной (см. _wilson_bubble_pressure_estimate) —
        # quick_bracket в _bisect_stability_boundary сам проверяет, годится
        # ли она, и откатывается на скан, если нет.
        if p_guess_bubble is None:
            p_guess_bubble = _wilson_bubble_pressure_estimate(composition, t_k)

        self._p_guess_bubble = p_guess_bubble
        self._p_guess_dew = p_guess_dew
        self.last_primary_status = NOT_PRESENT
        self.last_secondary_status = NOT_PRESENT
        self._stability_failed = False

    def _is_stable(self, p: float) -> bool | None:
        """
        Тест стабильности в точке `p`. Возвращает `True`/`False`, либо `None`,
        если тест не сошёлся за отведённое число итераций (`StopIterationError`
        из `TwoPhaseStabilityTest`, например, вблизи границы `p_max_bar`, если
        она выбрана слишком заниженной для этой температуры) — значение
        считается неопределённым, а не ошибкой, обрушивающей весь расчёт.
        """
        try:
            stab = TwoPhaseStabilityTest(self.composition, p, self.t_k)
            stab.calculate_phase_stability()
            return stab.stable
        except StopIterationError:
            self._stability_failed = True
            logger.warning(
                "Тест стабильности не сошёлся при P=%.4f бар (T=%.2f K)",
                p, self.t_k,
            )
            return None

    def _classify_point_type(self, p: float, nudge_sign: float) -> str:
        """
        Физическая классификация найденной точки насыщения — критерий
        Михельсена (тот же, что уже встроен в `TwoPhaseStabilityTest.
        k_vals_for_sat_point_calculation`): доминирует `S_v` (нарождается
        пар) → `'bubble'`, доминирует `S_l` (нарождается жидкость) → `'dew'`.

        Это НЕЗАВИСИМЫЙ от направления поиска шаг: раньше тип определялся
        только тем, сверху вниз или снизу вверх нашли точку, что физически
        неверно за критической точкой состава (см. докстринг модуля).

        Классифицируем не ровно в `p` (там флаг стабильности неоднозначен —
        точность самой бисекции `bracket_tol_rel`), а чуть сдвинувшись внутрь
        заведомо нестабильной (двухфазной) стороны — `nudge_sign` (+1 или -1)
        задаёт направление сдвига и известен вызывающему коду (для Bubble-
        поиска нестабильно ниже `p`, для Dew-поиска — выше).
        """
        p_test = p * (1.0 + nudge_sign * 10.0 * self.bracket_tol_rel)
        stab = TwoPhaseStabilityTest(self.composition, p_test, self.t_k)
        stab.calculate_phase_stability()
        return 'bubble' if (stab.S_v - 1.0) >= (stab.S_l - 1.0) else 'dew'

    def _bracket_is_tight(self, p_lo: float, p_hi: float) -> bool:
        """
        `True`, если ширина `[p_lo, p_hi]` уже не превышает допуск
        (относительный `bracket_tol_rel`, либо абсолютный `bracket_tol_abs` —
        что больше).

        `p_lo`/`p_hi` здесь означают "стабильный конец"/"нестабильный конец",
        а не "меньшее"/"большее" значение — для ветки Bubble (поиск сверху
        вниз) `p_lo` физически БОЛЬШЕ `p_hi`. Разность обязательно через
        `abs()` — без него для этой ветки результат всегда отрицателен и
        проверка ложно срабатывает на первой же итерации.
        """
        width = abs(p_hi - p_lo)
        scale = min(p_lo, p_hi)
        return width < max(self.bracket_tol_abs, self.bracket_tol_rel * scale)

    def _find_transition_scanning(
        self, p_from: float, p_to: float,
        initial_log_step: float = 0.05, max_iter: int = 40,
    ) -> tuple[float, float] | None:
        """
        Ищет ПЕРВУЮ (считая от `p_from`) смену флага стабильности, двигаясь
        от `p_from` к `p_to` геометрически растущим шагом по log(P): шаг
        начинается с `initial_log_step` (~5% по давлению) и удваивается на
        каждой итерации.

        Нужно именно сканирование, а не проверка только двух крайних точек:
        за критической точкой (не считаем её здесь — см. докстринг модуля)
        диапазон может быть "сэндвичем" стабильно-нестабильно-стабильно, и
        оба конца могут случайно совпасть по флагу, хотя переход внутри есть.

        Растущий (а не равномерный) шаг — сознательная оптимизация скорости:
        если перехода нет вовсе (типичный случай для Dew-ветки на бóльшей
        части диапазона T — второй точки насыщения физически нет), полный
        равномерный скан тратил фиксированное число проверок независимо от
        этого, а растущий шаг доходит до `p_to` всего за ~log2(диапазон)
        итераций. Одновременно у самой `p_from` разрешение остаётся мелким
        (~5%) — там, где по опыту и находится нижняя ретроградная точка (в
        первых процентах диапазона от p_min_bar). Пропустить переход можно
        только если реальная двухфазная область уже сузилась до ширины
        меньше текущего (растущего) шага — на практике это происходит
        только вплотную к самой критической точке, которую мы здесь и так
        не считаем (см. докстринг модуля).

        Возвращает `(p_stable, p_unstable)` — bracket вокруг найденного
        перехода (ближайшего к `p_from`) для последующей точной бисекции,
        либо `None`, если флаг ни разу не поменялся до `p_to` (перехода в
        этом диапазоне нет) или тест стабильности не сошёлся в самой `p_from`.
        """
        log_from = np.log(p_from)
        log_to = np.log(p_to)
        direction = 1.0 if log_to >= log_from else -1.0

        prev_p = p_from
        prev_flag = self._is_stable(prev_p)
        if prev_flag is None:
            return None

        step = initial_log_step
        log_p = log_from

        for _ in range(max_iter):
            log_p += direction * step
            reached_end = direction * (log_p - log_to) >= 0
            if reached_end:
                log_p = log_to

            p = float(np.exp(log_p))
            flag = self._is_stable(p)

            if flag is not None:
                if flag != prev_flag:
                    return (prev_p, p) if prev_flag else (p, prev_p)
                prev_p, prev_flag = p, flag

            if reached_end:
                break
            step *= 2.0

        return None

    def _bisect_stability_boundary(
        self,
        p_stable_end: float,
        p_unstable_end: float,
        quick_bracket: tuple[float, float] | None,
    ) -> float | None:
        """
        Бисекция по флагу стабильности между гарантированно стабильным
        (`p_stable_end`) и гарантированно нестабильным (`p_unstable_end`)
        концом полного диапазона. Если `quick_bracket` задан (узкая попытка
        вокруг guess с предыдущего шага по T) и на его концах флаги
        действительно разные — используется он (быстрее); иначе — безопасный
        откат на полный `[p_stable_end, p_unstable_end]`.

        Возвращает `None`, если даже на полном диапазоне флаг стабильности не
        меняется (перехода/второй точки насыщения в этом диапазоне нет) или
        тест стабильности не сошёлся на одном из контрольных концов.
        """
        p_lo = p_hi = None

        if quick_bracket is not None:
            p_a, p_b = quick_bracket
            if abs(p_a - p_stable_end) <= abs(p_b - p_stable_end):
                p_s, p_u = p_a, p_b
            else:
                p_s, p_u = p_b, p_a
            if self._is_stable(p_s) is True and self._is_stable(p_u) is False:
                p_lo, p_hi = p_s, p_u

        if p_lo is None:
            # За критической точкой (которую мы не считаем — см. докстринг
            # модуля) диапазон [p_stable_end, p_unstable_end] может содержать
            # НЕСКОЛЬКО переходов стабильности ("сэндвич" стабильно-нестабильно-
            # стабильно — типично для состава в ретроградной области: проверка
            # только двух крайних точек эту структуру не ловит вовсе, если оба
            # конца случайно оказываются стабильны). Сканируем от `p_stable_end`
            # в сторону `p_unstable_end`, находим БЛИЖАЙШИЙ к `p_stable_end`
            # переход — это либо искомая Bubble-точка (если сканируем сверху),
            # либо нижняя ретроградная Dew-точка (если сканируем снизу).
            scan_bracket = self._find_transition_scanning(p_stable_end, p_unstable_end)
            if scan_bracket is None:
                return None
            p_lo, p_hi = scan_bracket

        for _ in range(self.max_iter):
            if self._bracket_is_tight(p_lo, p_hi):
                break

            p_mid = (p_lo + p_hi) / 2.0
            flag_mid = self._is_stable(p_mid)

            if flag_mid is None:
                # Тест не сошёлся ровно в середине — отступаем ближе к уже
                # проверенной стабильной стороне и пробуем ещё раз один раз.
                p_mid = p_lo + (p_mid - p_lo) * 0.5
                flag_mid = self._is_stable(p_mid)
                if flag_mid is None:
                    logger.warning(
                        "Бисекция прервана: тест стабильности не сходится "
                        "вблизи P=%.4f бар (T=%.2f K)", p_mid, self.t_k,
                    )
                    # Возвращать середину всё ещё широкого bracket'а здесь
                    # опасно: это выглядело бы как валидная точка насыщения,
                    # хотя последняя проверка границы не состоялась.
                    return None

            if flag_mid:
                p_lo = p_mid
            else:
                p_hi = p_mid

        return (p_lo + p_hi) / 2.0

    def find_bubble_point(self) -> SaturationPoint | None:
        """
        Поиск первой (верхней) точки насыщения — переход стабильно(`p_max_bar`)
        → нестабильно(`p_min_bar`), т.е. сверху вниз. Несмотря на имя метода
        (сохранено по направлению поиска, как в `PhaseDiagram_v4.py`),
        физический тип возвращаемой точки может оказаться `'dew'` — см.
        `point_type` и докстринг модуля/`_classify_point_type`.
        """
        self._stability_failed = False
        self.last_primary_status = NOT_PRESENT
        quick_bracket = None
        if self._p_guess_bubble is not None:
            g = self._p_guess_bubble
            quick_bracket = (
                max(self.p_min_bar, g * 0.7),
                min(self.p_max_bar, g * 1.3),
            )

        result = self._bisect_stability_boundary(
            p_stable_end=self.p_max_bar,
            p_unstable_end=self.p_min_bar,
            quick_bracket=quick_bracket,
        )
        if result is None:
            if self._stability_failed:
                self.last_primary_status = STABILITY_NONCONVERGENCE
            logger.info("Верхняя точка насыщения не найдена (T=%.2f K)", self.t_k)
            return None

        try:
            point_type = self._classify_point_type(result, nudge_sign=-1.0)
        except StopIterationError:
            self.last_primary_status = STABILITY_NONCONVERGENCE
            logger.warning(
                "Классификация верхней точки не сошлась при P=%.4f бар (T=%.2f K)",
                result, self.t_k,
            )
            return None
        self.last_primary_status = OK
        logger.info(
            "Верхняя точка насыщения: P=%.4f бар (T=%.2f K), тип=%s",
            result, self.t_k, point_type,
        )
        return SaturationPoint(pressure=result, point_type=point_type)

    def find_dew_point(self, p_upper_bound: float) -> SaturationPoint | None:
        """
        Поиск второй (нижней) точки насыщения — переход стабильно(`p_min_bar`)
        → нестабильно(чуть ниже `p_upper_bound`), т.е. снизу вверх, в границах
        `[p_min_bar, p_upper_bound]` (`p_upper_bound` — обычно уже найденная на
        этой температуре первая/верхняя точка). Контрольная "нестабильная"
        точка берётся не строго в `p_upper_bound`, а с небольшим отступом
        внутрь — ровно на границе флаг стабильности неоднозначен из-за
        конечной точности, с которой была найдена сама `p_upper_bound`.

        Физически это почти всегда `'dew'` (нижняя ретроградная граница), но
        тип всё равно классифицируется тем же критерием, а не жёстко
        прописывается — см. `_classify_point_type`.
        """
        self._stability_failed = False
        self.last_secondary_status = NOT_PRESENT
        p_max = min(p_upper_bound, self.p_max_bar)
        p_unstable_ref = p_max - max(self.bracket_tol_abs * 10.0, p_max * self.bracket_tol_rel * 10.0)

        quick_bracket = None
        if self._p_guess_dew is not None:
            g = self._p_guess_dew
            quick_bracket = (
                max(self.p_min_bar, g * 0.7),
                min(p_unstable_ref, g * 1.3),
            )

        result = self._bisect_stability_boundary(
            p_stable_end=self.p_min_bar,
            p_unstable_end=p_unstable_ref,
            quick_bracket=quick_bracket,
        )
        if result is None:
            if self._stability_failed:
                self.last_secondary_status = STABILITY_NONCONVERGENCE
            logger.info("Нижняя (ретроградная) точка насыщения не найдена (T=%.2f K)", self.t_k)
            return None

        try:
            point_type = self._classify_point_type(result, nudge_sign=1.0)
        except StopIterationError:
            self.last_secondary_status = STABILITY_NONCONVERGENCE
            logger.warning(
                "Классификация нижней точки не сошлась при P=%.4f бар (T=%.2f K)",
                result, self.t_k,
            )
            return None
        self.last_secondary_status = OK
        if point_type != 'dew':
            logger.warning(
                "Нижняя точка насыщения P=%.4f бар (T=%.2f K) классифицирована "
                "как '%s' — физически неожиданно (ожидалась 'dew')",
                result, self.t_k, point_type,
            )
        else:
            logger.info(
                "Нижняя точка насыщения: P=%.4f бар (T=%.2f K), тип=%s",
                result, self.t_k, point_type,
            )
        return SaturationPoint(pressure=result, point_type=point_type)


def _calculate_envelope_chunk_worker(
    composition: Composition, t_c_chunk: list[float], p_max_bar: float, p_min_bar: float,
) -> list[tuple[float, SaturationPoint | None, SaturationPoint | None, str, str]]:
    """
    Считает НЕПРЕРЫВНЫЙ участок диапазона T последовательно, с continuation
    внутри участка (как `PhaseEnvelopeSSM.calculate()`) — вынесена на уровень
    модуля для корректной сериализации в joblib (тот же приём, что в
    `PhaseEnvelopeFromStability.py`/notebook `CriticalPointCalculatorParrallel`).

    "Холодный" поиск без guess'а нужен только для ПЕРВОЙ точки участка —
    остальные стартуют от результата предыдущей точки того же участка, как в
    обычном последовательном марше. Это специально сделано КУСКАМИ (а не по
    одной точке на задачу): раньше каждая точка была отдельной joblib-задачей
    без continuation вовсе, и почти каждая точка целиком платила "холодную"
    цену (полный скан) — крупные куски резко сокращают число таких холодных
    стартов (один на кусок, а не один на точку), а заодно и число пересылок
    `composition` между процессами (дороже для составов со многими
    компонентами — BIP-матрица nc×nc и т.п.). Участок получает одну
    независимую копию состава (`deep_copy=True`), чтобы избежать гонки за
    `composition.T` между воркерами.
    """
    comp = composition.new_composition(composition.composition, deep_copy=True)

    results = []
    prev_pb: float | None = None
    prev_pdew: float | None = None

    for t_c in t_c_chunk:
        t_k = float(t_c) + 273.15

        finder = SaturationPointSSM(
            comp, t_k, p_max_bar, p_min_bar,
            p_guess_bubble=prev_pb, p_guess_dew=prev_pdew,
        )
        primary = finder.find_bubble_point()
        p_upper_for_dew = primary.pressure if primary is not None else p_max_bar
        secondary = finder.find_dew_point(p_upper_for_dew)

        results.append((
            float(t_c), primary, secondary,
            finder.last_primary_status, finder.last_secondary_status,
        ))

        if primary is not None:
            prev_pb = primary.pressure
        if secondary is not None:
            prev_pdew = secondary.pressure

    return results


class PhaseEnvelopeSSM:
    """
    Фазовая огибающая методом последовательных приближений в диапазоне
    температур — порт `PhaseDiagram_v4.py::PhaseDiagram`, дополненный
    continuation по температуре (давление с предыдущего шага передаётся как
    стартовая попытка bracket'а для следующего — в `PhaseDiagram_v4.py` этого
    не было, каждый шаг стартовал заново от `p_max/2`). Есть и параллельная
    версия (`calculate_parallel`, через joblib) — диапазон T делится на
    непрерывные куски по числу воркеров, continuation работает внутри каждого
    куска, см. её докстринг за подробностями.
    """

    def __init__(
        self,
        composition: Composition,
        t_min_c: float,
        t_max_c: float,
        t_step_c: float,
        p_max_bar: float,
        p_min_bar: float = 0.1,
        dew_bubble_min_gap_rel: float = 0.01,
    ):
        self.composition = composition
        self.t_min_c = t_min_c
        self.t_max_c = t_max_c
        self.t_step_c = t_step_c
        self.p_max_bar = p_max_bar
        self.p_min_bar = p_min_bar
        # Ветка Dew ищется в границах [p_min_bar, Bubble] на той же T. Если на
        # этом составе/диапазоне T нет настоящей второй (ретроградной) точки
        # насыщения, единственный переход стабильности в этом интервале — это
        # сам Bubble, и ветка Dew может сойтись к почти тому же давлению
        # (отличие — доли процента, а не отдельная точка). Такие "почти-дубли"
        # отбрасываются (Dew -> NaN), если относительная разница меньше порога.
        self.dew_bubble_min_gap_rel = dew_bubble_min_gap_rel

        self.temps_c = np.arange(t_min_c, t_max_c + t_step_c / 2.0, t_step_c)
        self.results = {'Temp_C': [], 'Bubble_bar': [], 'Dew_upper_bar': [], 'Dew_lower_bar': []}
        self.diagnostics: list[EnvelopeDiagnostic] = []

    def _frame(self) -> pd.DataFrame:
        return attach_diagnostics(
            pd.DataFrame(self.results), self.diagnostics,
            method="ssm", fallback=None,
        )

    def _filter_near_duplicate_dew(
        self, t_c: float, primary: SaturationPoint | None, secondary: SaturationPoint | None,
    ) -> SaturationPoint | None:
        """Отбрасывает secondary (нижнюю точку), если она практически совпадает
        с primary (см. `dew_bubble_min_gap_rel`) — значит, второй точки
        насыщения на этой T физически нет, а найденный "дубль" — просто та же
        самая граница, подойдённая снизу."""
        if (
            primary is not None
            and secondary is not None
            and abs(secondary.pressure - primary.pressure) / primary.pressure < self.dew_bubble_min_gap_rel
        ):
            logger.info(
                "T=%.2f °C: нижняя точка (%.4f бар) слишком близко к верхней (%.4f бар) — "
                "считаем, что второй точки насыщения нет",
                t_c, secondary.pressure, primary.pressure,
            )
            return None
        return secondary

    def _append_result(
        self, t_c: float, primary: SaturationPoint | None, secondary: SaturationPoint | None,
    ) -> None:
        """
        Раскладывает найденные точки по трём физическим колонкам (см.
        докстринг класса) — по аналогии с `PhaseEnvelope.py`
        (`bubble_point`/`dew_point_upper`/`dew_point_lower`):
        `primary` (поиск сверху вниз) идёт в `Bubble_bar`, если физически это
        точка кипения, иначе (за критической точкой состава — см. докстринг
        модуля) в `Dew_upper_bar`; `secondary` (поиск снизу вверх, когда
        найден) — в `Dew_lower_bar`.
        """
        self.results['Temp_C'].append(float(t_c))

        if primary is None:
            self.results['Bubble_bar'].append(np.nan)
            self.results['Dew_upper_bar'].append(np.nan)
        elif primary.point_type == 'bubble':
            self.results['Bubble_bar'].append(primary.pressure)
            self.results['Dew_upper_bar'].append(np.nan)
        else:
            self.results['Bubble_bar'].append(np.nan)
            self.results['Dew_upper_bar'].append(primary.pressure)

        self.results['Dew_lower_bar'].append(secondary.pressure if secondary is not None else np.nan)

    def calculate(
        self,
        cancellation_token: CancellationToken | None = None,
        progress_callback: ProgressCallback | None = None,
    ) -> pd.DataFrame:
        """Последовательный марш по температуре с continuation (см. докстринг класса)."""
        logger.info(
            "Расчёт фазовой огибающей (SSM): T=[%.1f, %.1f] °C, шаг %.1f °C",
            self.t_min_c, self.t_max_c, self.t_step_c,
        )
        self.results = {'Temp_C': [], 'Bubble_bar': [], 'Dew_upper_bar': [], 'Dew_lower_bar': []}
        self.diagnostics = []

        prev_pb: float | None = None
        prev_pdew: float | None = None

        total = max(1, len(self.temps_c))
        for index, t_c in enumerate(self.temps_c, start=1):
            if cancellation_token is not None:
                cancellation_token.throw_if_cancelled()
            report_progress(progress_callback, (index - 1) / total,
                            f"Envelope temperature {index}/{total}")
            t_k = float(t_c) + 273.15

            finder = SaturationPointSSM(
                self.composition, t_k, self.p_max_bar, self.p_min_bar,
                p_guess_bubble=prev_pb, p_guess_dew=prev_pdew,
            )
            primary = finder.find_bubble_point()
            p_upper_for_dew = primary.pressure if primary is not None else self.p_max_bar
            secondary = finder.find_dew_point(p_upper_for_dew)
            secondary = self._filter_near_duplicate_dew(t_c, primary, secondary)
            secondary_status = finder.last_secondary_status if secondary is not None else (
                NOT_PRESENT
                if finder.last_secondary_status == OK
                else finder.last_secondary_status
            )

            self._append_result(t_c, primary, secondary)
            self.diagnostics.append({
                'Temp_C': float(t_c),
                'primary_status': finder.last_primary_status,
                'secondary_status': secondary_status,
            })

            if primary is not None:
                prev_pb = primary.pressure
            if secondary is not None:
                prev_pdew = secondary.pressure

            report_progress(progress_callback, index / total,
                            f"Envelope temperature {index}/{total}")

        return self._frame()

    def calculate_parallel(self, n_jobs: int = -1, backend: str = 'loky') -> pd.DataFrame:
        """
        То же самое, что `calculate()`, но диапазон T разбивается на
        непрерывные куски (по числу доступных воркеров), и каждый кусок
        считается параллельно — см. `_calculate_envelope_chunk_worker`.

        Внутри куска continuation работает как в обычном последовательном
        марше — "холодный" (без guess'а) поиск нужен только для первой точки
        КАЖДОГО КУСКА, а не для каждой точки диапазона (как было в первой
        версии параллельного расчёта, где каждая T была отдельной задачей без
        continuation вовсе — большие куски резко сокращают число холодных
        стартов и число пересылок `composition` между процессами). Это не
        затрагивает корректность: guess влияет только на скорость, safety-net
        (скан + бисекция по флагу) всё равно отрабатывает, если guess не
        подошёл — см. докстринг `_bisect_stability_boundary`.
        """
        n_workers = max(1, min(effective_n_jobs(n_jobs), len(self.temps_c)))
        chunks = [chunk for chunk in np.array_split(self.temps_c, n_workers) if len(chunk) > 0]

        logger.info(
            "Расчёт фазовой огибающей (SSM, параллельно, %d кусков на %d воркеров): "
            "T=[%.1f, %.1f] °C, шаг %.1f °C",
            len(chunks), n_workers, self.t_min_c, self.t_max_c, self.t_step_c,
        )

        chunk_results = Parallel(n_jobs=n_jobs, backend=backend)(
            delayed(_calculate_envelope_chunk_worker)(
                self.composition, chunk.tolist(), self.p_max_bar, self.p_min_bar,
            )
            for chunk in chunks
        )

        self.results = {'Temp_C': [], 'Bubble_bar': [], 'Dew_upper_bar': [], 'Dew_lower_bar': []}
        self.diagnostics = []
        for chunk_result in chunk_results:
            for t_c, primary, secondary, primary_status, secondary_status in chunk_result:
                secondary = self._filter_near_duplicate_dew(t_c, primary, secondary)
                if secondary is None and secondary_status == OK:
                    secondary_status = NOT_PRESENT
                self._append_result(t_c, primary, secondary)
                self.diagnostics.append({
                    'Temp_C': float(t_c),
                    'primary_status': primary_status,
                    'secondary_status': secondary_status,
                })

        return self._frame()

    def get_data(self) -> pd.DataFrame:
        return self._frame()

    def plot(self, show: bool = True, save_path: str | None = None):
        """
        P-T диаграмма фазовой огибающей как ОДНОЙ непрерывной линии (без
        маркеров точек).

        `Bubble_bar`/`Dew_upper_bar` — взаимодополняющие половины одной и той
        же физической кривой (см. докстринг класса: за критической точкой
        состава поиск "сверху вниз" находит уже не bubble, а верхнюю точку
        росы) — объединяются в одну верхнюю ветку через `combine_first`, без
        разрыва в точке перехода.

        `Dew_lower_bar` (нижняя ретроградная точка) физически смыкается с
        верхней веткой у критической точки состава — пристыковывается к ней в
        обратном порядке по температуре, замыкая контур одной линией вместо
        трёх стилистически разных кусков. Из-за отбрасывания "почти-дублей"
        (`dew_bubble_min_gap_rel`, см. `_filter_near_duplicate_dew`) сама
        точка смыкания в данных отсутствует — линия соединяет ближайшие
        оставшиеся точки напрямую, без дополнительного сглаживания.
        """
        df = pd.DataFrame(self.results).sort_values('Temp_C').reset_index(drop=True)
        if df.empty:
            logger.warning("Нет данных для построения — вызовите calculate() перед plot().")
            return None

        fig, ax = plt.subplots(figsize=(10, 6))

        upper = df['Bubble_bar'].combine_first(df['Dew_upper_bar'])
        upper_mask = upper.notna()
        upper_t = df.loc[upper_mask, 'Temp_C'].to_numpy()
        upper_p = upper[upper_mask].to_numpy()

        lower = df.dropna(subset=['Dew_lower_bar']).sort_values('Temp_C', ascending=False)
        lower_t = lower['Temp_C'].to_numpy()
        lower_p = lower['Dew_lower_bar'].to_numpy()

        envelope_t = np.concatenate([upper_t, lower_t])
        envelope_p = np.concatenate([upper_p, lower_p])

        if envelope_t.size:
            ax.plot(envelope_t, envelope_p, 'b-', linewidth=2)

        ax.set_xlabel('Temperature, °C')
        ax.set_ylabel('Pressure, bar')
        ax.set_title('Phase envelope')
        ax.grid(True, linestyle='--', alpha=0.6)
        fig.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info("График сохранён: %s", save_path)
        if show:
            plt.show()

        return fig
