"""
Фазовая огибающая методом Ньютона по аналитическим производным летучести.

Второй (независимый) способ построить огибающую в дополнение к
`PhaseEnvelopeSuccessiveSubstitution.py` (SSM) — тот модуль НЕ меняется вообще,
остаётся инструментом для критической/ретроградной области (см. ниже).

Мотивация: `BubblePointPressure.py::BubblePointCalculator` и
`DewPressure.py::DewPointCalculator` уже решают точку насыщения методом Ньютона
по `BrusilovskiyEOS.calc_d_log_phi_i_dp` и на практике совпадают с PVTSim
почти точно (на составе KRSNL при T=110°C: 185.05 бар против 185.04 у PVTSim,
у SSM в той же точке — заметно грубее). Здесь эти два калькулятора
используются как "движок точки" внутри march'а по температуре (архитектурно
— тот же приём, что у `PhaseEnvelopeSSM`: continuation P между соседними T,
primary/secondary точка на T, тот же формат результата), только вместо
бисекции по флагу стабильности — Ньютон с continuation-guess.

Производной `d ln(phi)/dT` в проекте нет (только `calc_d_log_phi_i_dp` — по
P), поэтому буквальный двумерный алгоритм построения огибающей по спецификации
переменной (как в документации PVTsim: совместный Ньютон по T и P с
экстраполяцией по производным) здесь не реализован — вместо этого T задаётся
явной сеткой (как в SSM), а Ньютон на каждом шаге ищет только P.

ВАЖНОЕ ОГРАНИЧЕНИЕ (обнаружено эмпирически, определяет весь дизайн модуля):
`BubblePointCalculator` (уравнение "bubble": Σzᵢ·Kᵢ−1=0) и
`DewPointCalculator` (уравнение "dew": Σzᵢ/Kᵢ−1=0) — РАЗНЫЕ уравнения,
совпадающие только в самой критической точке состава. За критической точкой
bubble-уравнение всё равно может СОЙТИСЬ (с высокой точностью по своему
собственному F) к физически невернóму корню — на составе KRSNL за критической
точкой (~382°C) numeric-эксперимент показал 87.6 бар вместо истинных 126.8
(по SSM) при T=440°C, и калькулятор при этом рапортует `converged=True`.
Поэтому каждая найденная методом Ньютона точка здесь ОБЯЗАТЕЛЬНО проверяется
независимо — `TwoPhaseStabilityTest` по ОБЕ стороны от найденного `p` (см.
`SaturationPointNewton._verify`): принимается только если по одну сторону
система стабильна, а по другую — нестабильна (это и есть подтверждение, что
`p` — граница фаз), и физическая классификация (S_v vs S_l) совпадает с типом
уравнения, которым точка была найдена. Иначе — корень отбраковывается,
пробуется другой тип уравнения.

Проверка ДОЛЖНА быть именно двусторонней, а не только "вглубь двухфазной
области" (как в `SaturationPointSSM._classify_point_type` — для SSM
достаточно одной стороны, там `p` уже гарантированно граница, найденная
бисекцией; здесь же `p` — сырой корень уравнения Ньютона, который ещё не
факт, что вообще граница). Обнаружено эмпирически (состав KRSNL, T=390°C,
dew-уравнение): даже стартуя Ньютон ПРЯМО с истинного корня (~169.6 бар,
подтверждён SSM и прямым сканированием стабильности), первая же итерация
уводит его к P≈161.2 бар — точке, которая ВНУТРИ двухфазной области (не
на границе): при одностороннем сдвиге вниз это выглядело как "нестабильно,
похоже на dew" и проходило старую (одностороннюю) проверку, но сдвиг вверх
от 161.2 ТОЖЕ показывает нестабильно — то есть 161.2 не граница вообще.
Причина — судя по всему, в самом `DewPointCalculator`/`BubblePointCalculator`:
их аналитическая производная `dF/dP` почти на каждой итерации получается
"неправильного" знака и принудительно инвертируется веткой "Sign correction"
(см. побочную находку в конце этого докстринга) — вблизи критической точки
это может уводить шаг не туда даже от идеального старта. Это свойство
существующих калькуляторов (не трогаются в рамках этой задачи), но именно
поэтому проверка здесь обязана быть двусторонней и независимой, а не
доверять `converged=True`.

С этой (исправленной) проверкой march уверенно и быстро идёт от низких T
через bubble-ветку, корректно переключается на dew (upper) ровно там же, где
и SSM находит переход, и лишь у самой критической точки состава перестаёт
находить верифицируемый корень — что ожидаемо: сам SSM-модуль уже
документирует, что "метод Ньютона... теряет устойчивость вблизи критической
точки". Поэтому:

- march ОСТАНАВЛИВАЕТСЯ (не пытается героически восстановиться холодным
  стартом — численно проверено, что Wilson-догадка для dew-upper за критикой
  на этом составе сходится к спонтанному ~500-600 бар, а не к истинным
  ~100-170) при первом неверифицированном шаге, заполняя остаток запрошенного
  диапазона T NaN.

Нижняя ретроградная dew-точка (`Dew_lower_bar`) устроена иначе: у неё нет
надёжного холодного старта по Вильсону (проверено — сходится не к тому
корню), НО, в отличие от верхней ветки у самой критической точки, вдали от
неё это хорошо обусловленное Ньютоновское уравнение — если дать ему ХОТЯ БЫ
ОДНУ хорошую стартовую точку, дальше по T оно само бежит быстро и устойчиво
(численно проверено: от разового seed при T=400°C(3.17 бар) до T=465°C(21.5
бар) на составе KRSNL — 14 точек за 0.67 сек, 3-21 итераций на точку,
совпадение с SSM до 4-5 значащих цифр). Поэтому `calculate()` (только
последовательный вариант, не `calculate_parallel()`) после остановки
основного march'а ОДИН РАЗ прибегает к `SaturationPointSSM.find_dew_point`
(bootstrap, `_bootstrap_lower_dew_seed`) — пробует несколько соседних T сразу
после точки остановки (ретроградная область только-только открывается там,
поэтому первые попытки могут быть слишком узкими и не дать разрешимого
результата), находит один seed, и дальше продолжает уже чистым Ньютоном
(`SaturationPointNewton.find_dew_point`) до конца запрошенного диапазона —
с той же верификацией и тем же бюджетом времени на шаг, что и у основной
ветки. Это единственное место в модуле, где `SaturationPointSSM` вообще
используется (сам SSM-модуль по-прежнему не редактируется ни строкой) — как
редкий, ограниченный по числу вызовов "затравочный" инструмент, а не как
основной путь расчёта.

Для верхней ветки (`_bootstrap_and_march_upper_dew`) СИММЕТРИЧНАЯ идея —
затравка через `SaturationPointSSM.find_bubble_point`, потом Ньютон — тоже
реализована, но, в отличие от нижней ветки, не гарантирует продвижения СРАЗУ
ЗА критической точкой: это не вопрос качества затравки (двусторонняя
верификация теперь надёжно отбраковывает плохие корни — см. выше), а
свойство самих `BubblePointCalculator`/`DewPointCalculator` там же — их
Ньютон-шаг может не сойтись к границе вообще ни от какого guess'а вплотную к
критике (та же причина, что и в находке про "Sign correction" выше). Поэтому
`max_cycles=1` (не пытаемся героически перезатравливаться много раз подряд
именно в этой полосе — предыдущая версия с `max_cycles=5` пробовала и просто
теряла время, скатываясь по стоимости к самой SSM, не получая взамен ничего
воспроизводимо верного). Результат: `Dew_upper_bar` сразу за критической
точкой состава честно остаётся NaN шириной в несколько шагов T (столько же,
сколько занимает основной march), пока верхняя ветка снова не станет хорошо
обусловленной — эта полоса покрывается только `PhaseEnvelopeSSM`.

Вблизи критической точки отдельные шаги могут резко замедляться (до
секунды и больше на точку) — не из-за числа итераций Ньютона, а из-за
внутренней логики `BubblePointCalculator._find_unstable_region`/
`_escape_trivial_solution` (сама делает вложенную бисекцию по
`TwoPhaseStabilityTest`, если continuation-guess ошибочно классифицируется
как "стабильный" — вблизи критики двухфазная область узкая, и это происходит
чаще). Не баг и не трогается (калькуляторы не меняются) — но чтобы march не
"залипал" на подходе к критике на много секунд, здесь добавлен бюджет времени
на шаг: несколько подряд медленных шагов (дольше `slow_point_budget_s`)
останавливают march проактивно, не дожидаясь первого явного отказа.

`tol` калькуляторов по умолчанию ослаблен до `1e-4` (их собственные дефолты
`1e-8`/`1e-10` избыточно точны — практическое расхождение с PVTSim и так
~0.03%, см. докстринг SSM) — это на порядок сокращает число итераций Ньютона
в "спокойной" зоне (с ~40-70 до ~10-20) без потери практической точности.

ПОБОЧНАЯ НАХОДКА (не часть этой задачи, объясняет ограничение выше): у
`BubblePointCalculator`/`DewPointCalculator` производная `dF/dP` почти на
каждой итерации получается "не того" знака и принудительно инвертируется
веткой `if dF_dP > 0: dF_dP = -abs(dF_dP)` — судя по названию ветки
("Sign correction"), задумывалась как редкий edge case, а по факту (по логам)
срабатывает почти всегда, даже вдали от критической точки, где для этого нет
явной физической причины. Вдали от критики это не портит результат
(подтверждено многократно: см. T=110°C выше) — просто больше итераций, чем
у "настоящего" Ньютона. Но именно вплотную к критической точке состава эта
принудительная инверсия, судя по всему, может увести шаг не туда даже от
идеального guess'а (см. пример с T=390°C выше) — отсюда и ограничение на
`Dew_upper_bar` сразу за критикой. Не исправляю в рамках этой задачи (файлы
не мои для правки без отдельного запроса).
"""

import logging
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from joblib import Parallel, delayed, effective_n_jobs

from calc_core.Composition.Composition import Composition
from calc_core.PhaseDiagram.BubblePointPressure import BubblePointCalculator
from calc_core.PhaseDiagram.DewPressure import DewPointCalculator
from calc_core.PhaseDiagram.PhaseEnvelopeSuccessiveSubstitution import SaturationPoint, SaturationPointSSM
from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest

logger = logging.getLogger(__name__)


class SaturationPointNewton:
    """
    Поиск давления насыщения на фиксированной температуре методом Ньютона по
    аналитическим производным летучести — движок точки: `BubblePointCalculator`
    (уравнение "bubble") / `DewPointCalculator` (уравнение "dew"), оба без
    изменений. Каждый найденный корень независимо верифицируется через
    `TwoPhaseStabilityTest` (см. докстринг модуля) — без этого Ньютон может
    тихо сойтись к физически неверному корню за критической точкой состава.
    """

    def __init__(
        self,
        composition: Composition,
        t_k: float,
        p_max_bar: float,
        p_min_bar: float = 0.1,
        p_guess_bubble: float | None = None,
        p_guess_dew: float | None = None,
        branch_hint: str = 'bubble',
        tol: float = 1e-4,
        max_iter: int = 60,
        verify_nudge_rel: float = 1e-3,
    ):
        self.composition = composition
        self.t_k = t_k
        self.p_max_bar = p_max_bar
        self.p_min_bar = p_min_bar
        self.p_guess_bubble = p_guess_bubble
        self.p_guess_dew = p_guess_dew
        self.branch_hint = branch_hint if branch_hint in ('bubble', 'dew') else 'bubble'
        self.tol = tol
        self.max_iter = max_iter
        self.verify_nudge_rel = verify_nudge_rel

    def _solve_bubble(self, guess: float | None) -> tuple[float | None, bool]:
        calc = BubblePointCalculator(
            self.composition, self.t_k, P_guess=guess, tol=self.tol, max_iter=self.max_iter,
        )
        p = calc.calculate()
        return p, calc.converged

    def _solve_dew(self, guess: float | None, dew_point_type: str) -> tuple[float | None, bool]:
        calc = DewPointCalculator(
            self.composition, self.t_k, dew_point_type=dew_point_type,
            P_guess=guess, tol=self.tol, max_iter=self.max_iter,
        )
        p = calc.calculate()
        return p, calc.converged

    def _verify(self, p: float, nudge_sign: float) -> str | None:
        """
        Независимая физическая проверка корня, найденного Ньютоном —
        ДВУСТОРОННИЙ сдвиг от `p` и `TwoPhaseStabilityTest` по обе стороны.

        Изначально здесь была только ОДНОСТОРОННЯЯ проверка (сдвиг только
        вглубь предполагаемой двухфазной области, `nudge_sign`) — тот же
        приём, что `SaturationPointSSM._classify_point_type`. Для SSM этого
        достаточно: там `p` уже гарантированно лежит НА границе (найдена
        бисекцией по флагу стабильности), односторонний сдвиг там нужен только
        для классификации bubble/dew. Здесь же `p` — сырой корень уравнения
        Ньютона, который ещё предстоит подтвердить как границу вообще, а не
        просто точку ВНУТРИ двухфазной области — а уравнение (Σzᵢ·Kᵢ−1=0 или
        Σzᵢ/Kᵢ−1=0) может иметь и такие "внутренние" корни. Обнаружено
        эмпирически: на составе KRSNL при T=390°C Ньютон (dew-уравнение)
        сходится к P=161.17 бар, что ПРОХОДИЛО одностороннюю проверку (сдвиг
        вниз — нестабильно), но при сдвиге ВВЕРХ тоже нестабильно (истинная
        граница у этого состава при этой T — P≈169.5 бар) — то есть 161.17
        просто лежит внутри двухфазной области, а не на её границе.
        Одностороння проверка такой случай не ловит вообще.

        Returns
        -------
        str | None
            `'bubble'`/`'dew'` — физическая классификация точки, если по одну
            сторону от `p` система нестабильна, а по другую — стабильна (это
            и есть подтверждение, что `p` — граница фаз, а не случайный
            внутренний корень); `None`, если обе стороны нестабильны (корень
            внутри двухфазной области — не граница) или обе стабильны (сдвиг
            вообще не задел двухфазную область).
        """
        p_inside = p * (1.0 + nudge_sign * 10.0 * self.verify_nudge_rel)
        p_outside = p * (1.0 - nudge_sign * 10.0 * self.verify_nudge_rel)

        stab_inside = TwoPhaseStabilityTest(self.composition, p_inside, self.t_k)
        stab_inside.calculate_phase_stability()
        if stab_inside.stable:
            return None

        stab_outside = TwoPhaseStabilityTest(self.composition, p_outside, self.t_k)
        stab_outside.calculate_phase_stability()
        if not stab_outside.stable:
            return None

        return 'bubble' if (stab_inside.S_v - 1.0) >= (stab_inside.S_l - 1.0) else 'dew'

    def find_bubble_point(self) -> SaturationPoint | None:
        """
        Поиск первой (верхней) точки насыщения. Имя метода сохранено по
        аналогии с `SaturationPointSSM.find_bubble_point` — физический тип
        результата может оказаться `'dew'` (верхняя ретроградная точка за
        критической точкой состава).

        Пробует уравнение типа `branch_hint` первым (обычно — тип, найденный
        на предыдущем шаге T при continuation), при неудаче/провале
        верификации — второй тип. `None`, если ни один тип не даёт
        верифицируемого корня.
        """
        order = ['bubble', 'dew'] if self.branch_hint == 'bubble' else ['dew', 'bubble']

        for kind in order:
            if kind == 'bubble':
                p, converged = self._solve_bubble(self.p_guess_bubble)
            else:
                p, converged = self._solve_dew(self.p_guess_dew, 'upper')

            if not converged or p is None:
                continue
            if not (self.p_min_bar < p < self.p_max_bar):
                continue

            label = self._verify(p, nudge_sign=-1.0)
            if label != kind:
                logger.debug(
                    "Newton: корень P=%.4f бар (T=%.2f K) по '%s'-уравнению не "
                    "прошёл верификацию (метка=%s) — отбраковка",
                    p, self.t_k, kind, label,
                )
                continue

            logger.info(
                "Newton: верхняя точка насыщения P=%.4f бар (T=%.2f K), тип=%s",
                p, self.t_k, label,
            )
            return SaturationPoint(pressure=float(p), point_type=label)

        logger.info("Newton: верхняя точка насыщения не найдена/не верифицирована (T=%.2f K)", self.t_k)
        return None

    def find_dew_point(self, p_upper_bound: float) -> SaturationPoint | None:
        """
        Поиск второй (нижней, ретроградной) точки насыщения — только
        dew-уравнение (`dew_point_type='lower'`), только по continuation-guess
        (`p_guess_dew`). Без явного guess практически не сходится к верному
        корню (см. докстринг модуля) — в отличие от `SaturationPointSSM`, тут
        нет бисекции/скана, которые могли бы обнаружить эту точку "вслепую".
        """
        if self.p_guess_dew is None:
            logger.debug(
                "Newton: нижняя точка насыщения не ищется без guess (T=%.2f K) — "
                "у метода нет надёжного холодного старта для этой ветки",
                self.t_k,
            )
            return None

        p, converged = self._solve_dew(self.p_guess_dew, 'lower')
        if not converged or p is None:
            return None
        if not (self.p_min_bar < p < min(p_upper_bound, self.p_max_bar)):
            return None

        label = self._verify(p, nudge_sign=1.0)
        if label != 'dew':
            logger.debug(
                "Newton: нижняя точка P=%.4f бар (T=%.2f K) не прошла "
                "верификацию (метка=%s) — отбраковка", p, self.t_k, label,
            )
            return None

        logger.info("Newton: нижняя точка насыщения P=%.4f бар (T=%.2f K)", p, self.t_k)
        return SaturationPoint(pressure=float(p), point_type='dew')


def _calculate_envelope_chunk_worker(
    composition: Composition, t_c_chunk: list[float], p_max_bar: float, p_min_bar: float,
    tol: float, max_iter: int, slow_point_budget_s: float, max_consecutive_slow_points: int,
) -> list[tuple[float, SaturationPoint | None, SaturationPoint | None]]:
    """
    Считает НЕПРЕРЫВНЫЙ участок диапазона T последовательно, с continuation
    внутри участка — вынесена на уровень модуля для сериализации в joblib
    (тот же приём, что в `PhaseEnvelopeSSM._calculate_envelope_chunk_worker`).

    Холодный старт (Wilson, через `p_guess_bubble=None`) нужен только для
    первой точки участка. Если она уже лежит за критической/ретроградной
    областью состава, холодный старт её, скорее всего, не найдёт (см.
    докстринг модуля) — тогда весь участок вернётся пустым.
    """
    comp = composition.new_composition(composition.composition, deep_copy=True)

    results = []
    prev_pb: float | None = None
    prev_pdew: float | None = None
    prev_type = 'bubble'
    slow_streak = 0

    for t_c in t_c_chunk:
        t_k = float(t_c) + 273.15
        t0 = time.time()

        finder = SaturationPointNewton(
            comp, t_k, p_max_bar, p_min_bar,
            p_guess_bubble=prev_pb, p_guess_dew=prev_pdew,
            branch_hint=prev_type, tol=tol, max_iter=max_iter,
        )
        primary = finder.find_bubble_point()
        p_upper_for_dew = primary.pressure if primary is not None else p_max_bar
        secondary = finder.find_dew_point(p_upper_for_dew)

        elapsed = time.time() - t0
        slow_streak = slow_streak + 1 if elapsed > slow_point_budget_s else 0

        if primary is None or slow_streak >= max_consecutive_slow_points:
            results.append((float(t_c), None, None))
            break

        results.append((float(t_c), primary, secondary))
        prev_pb, prev_type = primary.pressure, primary.point_type
        if secondary is not None:
            prev_pdew = secondary.pressure

    return results


class PhaseEnvelopeNewton:
    """
    Фазовая огибающая методом Ньютона по производным летучести — второй,
    независимый от `PhaseEnvelopeSSM` способ построить огибающую (см.
    докстринг модуля за мотивацией/ограничениями). Тот же интерфейс и формат
    результата, что у `PhaseEnvelopeSSM` (для прямого сравнения/наложения на
    графике). Основная (bubble/dew upper) ветка ЧИСТО ОСТАНАВЛИВАЕТСЯ у
    критической точки состава, вместо того чтобы пытаться пройти её (как это
    делает SSM) — но `calculate()` (не `calculate_parallel()`) следом
    достраивает нижнюю ретроградную ветку (`Dew_lower_bar`) отдельным
    проходом: одна затравочная точка через `SaturationPointSSM`, дальше —
    снова чистый Ньютон (см. `_bootstrap_and_march_lower_dew`).
    """

    def __init__(
        self,
        composition: Composition,
        t_min_c: float,
        t_max_c: float,
        t_step_c: float,
        p_max_bar: float,
        p_min_bar: float = 0.1,
        tol: float = 1e-4,
        max_iter: int = 60,
        dew_bubble_min_gap_rel: float = 0.01,
        slow_point_budget_s: float = 0.5,
        max_consecutive_slow_points: int = 2,
        bootstrap_lower_dew: bool = True,
        lower_dew_bootstrap_attempts: int = 6,
        bootstrap_bracket_tol_rel: float = 0.03,
    ):
        self.composition = composition
        self.t_min_c = t_min_c
        self.t_max_c = t_max_c
        self.t_step_c = t_step_c
        self.p_max_bar = p_max_bar
        self.p_min_bar = p_min_bar
        self.tol = tol
        self.max_iter = max_iter
        self.dew_bubble_min_gap_rel = dew_bubble_min_gap_rel
        self.slow_point_budget_s = slow_point_budget_s
        self.max_consecutive_slow_points = max_consecutive_slow_points
        # Затравка нижней ретроградной ветки через SSM (см. докстринг модуля) —
        # запускается только в calculate() (не в calculate_parallel()), только
        # один раз за march, и только если основная ветка вообще покинула
        # bubble-область в запрошенном диапазоне T.
        self.bootstrap_lower_dew = bootstrap_lower_dew
        self.lower_dew_bootstrap_attempts = lower_dew_bootstrap_attempts
        # Затравочные вызовы SSM используют НАМНОГО более грубый bracket_tol_rel,
        # чем дефолт SSM (1e-3) — нам нужна только стартовая точка для Ньютона,
        # который сам уточнит её на порядки точнее; на этом составе вблизи
        # критической точки это отличие — разница между ~0.6 сек и ~5+ сек на
        # один затравочный вызов (см. докстринг модуля).
        self.bootstrap_bracket_tol_rel = bootstrap_bracket_tol_rel

        self.temps_c = np.arange(t_min_c, t_max_c + t_step_c / 2.0, t_step_c)
        self.results = {'Temp_C': [], 'Bubble_bar': [], 'Dew_upper_bar': [], 'Dew_lower_bar': []}

    def _filter_near_duplicate_dew(
        self, t_c: float, primary: SaturationPoint | None, secondary: SaturationPoint | None,
    ) -> SaturationPoint | None:
        """Отбрасывает secondary, если она практически совпадает с primary —
        см. `PhaseEnvelopeSSM._filter_near_duplicate_dew` (та же логика,
        продублирована здесь, чтобы не импортировать приватную зависимость
        из SSM-модуля)."""
        if (
            primary is not None
            and secondary is not None
            and abs(secondary.pressure - primary.pressure) / primary.pressure < self.dew_bubble_min_gap_rel
        ):
            return None
        return secondary

    def _append_result(
        self, t_c: float, primary: SaturationPoint | None, secondary: SaturationPoint | None,
    ) -> None:
        """Раскладывает найденные точки по трём физическим колонкам — см.
        `PhaseEnvelopeSSM._append_result` (та же логика)."""
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

    def _append_nan_rows(self, temps_c) -> None:
        for t_c in temps_c:
            self.results['Temp_C'].append(float(t_c))
            self.results['Bubble_bar'].append(np.nan)
            self.results['Dew_upper_bar'].append(np.nan)
            self.results['Dew_lower_bar'].append(np.nan)

    def calculate(self) -> pd.DataFrame:
        """
        Последовательный марш по температуре с continuation. Останавливается
        (заполняя остаток диапазона NaN) при первом неверифицированном шаге
        или после нескольких подряд аномально медленных шагов (см. докстринг
        модуля) — обычно это означает подход к критической точке состава, для
        которой предназначен `PhaseEnvelopeSSM`.
        """
        logger.info(
            "Расчёт фазовой огибающей (Newton): T=[%.1f, %.1f] °C, шаг %.1f °C",
            self.t_min_c, self.t_max_c, self.t_step_c,
        )
        self.results = {'Temp_C': [], 'Bubble_bar': [], 'Dew_upper_bar': [], 'Dew_lower_bar': []}

        prev_pb: float | None = None
        prev_pdew: float | None = None
        prev_type = 'bubble'
        slow_streak = 0

        for i, t_c in enumerate(self.temps_c):
            t_k = float(t_c) + 273.15
            t0 = time.time()

            finder = SaturationPointNewton(
                self.composition, t_k, self.p_max_bar, self.p_min_bar,
                p_guess_bubble=prev_pb, p_guess_dew=prev_pdew,
                branch_hint=prev_type, tol=self.tol, max_iter=self.max_iter,
            )
            primary = finder.find_bubble_point()
            p_upper_for_dew = primary.pressure if primary is not None else self.p_max_bar
            secondary = finder.find_dew_point(p_upper_for_dew)
            secondary = self._filter_near_duplicate_dew(t_c, primary, secondary)

            elapsed = time.time() - t0
            slow_streak = slow_streak + 1 if elapsed > self.slow_point_budget_s else 0

            if primary is None:
                logger.info(
                    "Newton march остановлен на T=%.2f °C (точка не верифицирована) — "
                    "вероятно, окрестность критической точки состава; для этой зоны "
                    "используйте PhaseEnvelopeSSM", t_c,
                )
                self._append_nan_rows(self.temps_c[i:])
                break

            if slow_streak >= self.max_consecutive_slow_points:
                logger.info(
                    "Newton march остановлен на T=%.2f °C (%d шагов подряд дольше "
                    "%.1f сек) — вероятно, окрестность критической точки состава; "
                    "для этой зоны используйте PhaseEnvelopeSSM",
                    t_c, slow_streak, self.slow_point_budget_s,
                )
                self._append_result(t_c, primary, secondary)
                self._append_nan_rows(self.temps_c[i + 1:])
                break

            self._append_result(t_c, primary, secondary)
            prev_pb, prev_type = primary.pressure, primary.point_type
            if secondary is not None:
                prev_pdew = secondary.pressure

        if self.bootstrap_lower_dew:
            # Сначала верхняя ветка (Bubble/Dew_upper) — нижняя ветка сверяется
            # с ней при отбраковке дублей (см. ниже), поэтому порядок важен.
            self._bootstrap_and_march_upper_dew()
            self._bootstrap_and_march_lower_dew()

        return pd.DataFrame(self.results)

    def _bootstrap_and_march_upper_dew(self, max_cycles: int = 1) -> None:
        """
        Достраивает верхнюю ветку (`Bubble_bar`/`Dew_upper_bar`) за точкой, где
        основной march остановился — тем же приёмом, что и нижняя ретроградная
        ветка (см. `_bootstrap_and_march_lower_dew`): затравка через
        `SaturationPointSSM.find_bubble_point` (несмотря на имя — как и в самом
        SSM, физический тип найденной точки может оказаться `'dew'` — обычный
        случай здесь, раз мы уже за критической точкой состава), потом снова
        чистый Ньютон.

        В отличие от нижней ветки, здесь ОДНОЙ затравки обычно недостаточно:
        основной march останавливается именно там, где считать труднее всего
        (упор в саму критическую точку — Z_L≈Z_V), и "трудная зона" вокруг неё
        имеет заметную ширину по T — Ньютон после затравки типично проходит
        несколько точек и снова упирается в тот же эффект чуть дальше. Поэтому
        весь процесс — ЦИКЛ (до `max_cycles` раз): затравка SSM → continuation
        Ньютоном, пока не остановится → если дошли до конца диапазона, стоп;
        если застряли раньше — затравка SSM ещё раз с места остановки, и так
        далее, пока либо не дойдём до конца, либо очередная затравка не
        продвинет вперёд ни на шаг (тупик — дальше не пытаемся).
        """
        bubble = self.results['Bubble_bar']
        dew_upper = self.results['Dew_upper_bar']

        for _cycle in range(max_cycles):
            start_idx = next(
                (i for i, (b, d) in enumerate(zip(bubble, dew_upper)) if pd.isna(b) and pd.isna(d)),
                None,
            )
            if start_idx is None:
                return  # весь диапазон уже покрыт

            seed_idx, seed_p, seed_type = None, None, None
            attempts = min(self.lower_dew_bootstrap_attempts, len(self.temps_c) - start_idx)
            for j in range(start_idx, start_idx + attempts):
                t_c = self.temps_c[j]
                t_k = float(t_c) + 273.15
                ssm_finder = SaturationPointSSM(
                    self.composition, t_k, self.p_max_bar, self.p_min_bar,
                    bracket_tol_rel=self.bootstrap_bracket_tol_rel,
                )
                pt = ssm_finder.find_bubble_point()
                if pt is not None:
                    seed_idx, seed_p, seed_type = j, pt.pressure, pt.point_type
                    logger.info(
                        "Newton: верхняя ветка затравлена через SSM на T=%.2f °C, "
                        "P=%.4f бар, тип=%s", t_c, seed_p, seed_type,
                    )
                    break

            if seed_idx is None:
                logger.info(
                    "Newton: не удалось затравить верхнюю ветку через SSM "
                    "в первых %d точках после T=%.2f °C — останавливаюсь",
                    attempts, self.temps_c[start_idx],
                )
                return

            prev_p, prev_type = seed_p, seed_type
            slow_streak = 0
            progressed = False
            for j in range(seed_idx, len(self.temps_c)):
                t_c = self.temps_c[j]
                t_k = float(t_c) + 273.15
                t0 = time.time()

                finder = SaturationPointNewton(
                    self.composition, t_k, self.p_max_bar, self.p_min_bar,
                    p_guess_bubble=prev_p if prev_type == 'bubble' else None,
                    p_guess_dew=prev_p if prev_type == 'dew' else None,
                    branch_hint=prev_type, tol=self.tol, max_iter=self.max_iter,
                )
                pt = finder.find_bubble_point()

                elapsed = time.time() - t0
                slow_streak = slow_streak + 1 if elapsed > self.slow_point_budget_s else 0

                if pt is None:
                    logger.info(
                        "Newton: верхняя ветка (после затравки) остановлена на "
                        "T=%.2f °C (точка не верифицирована)", t_c,
                    )
                    break

                if pt.point_type == 'bubble':
                    bubble[j] = pt.pressure
                else:
                    dew_upper[j] = pt.pressure
                prev_p, prev_type = pt.pressure, pt.point_type
                progressed = True

                if slow_streak >= self.max_consecutive_slow_points:
                    logger.info(
                        "Newton: верхняя ветка (после затравки) остановлена на "
                        "T=%.2f °C (%d шагов подряд дольше %.1f сек)",
                        t_c, slow_streak, self.slow_point_budget_s,
                    )
                    break
            else:
                return  # дошли до конца диапазона без остановки — готово

            if not progressed:
                logger.info(
                    "Newton: затравка на T=%.2f °C не продвинула верхнюю ветку "
                    "ни на шаг — останавливаюсь", self.temps_c[seed_idx],
                )
                return

    def _bootstrap_and_march_lower_dew(self) -> None:
        """
        Достраивает нижнюю ретроградную dew-ветку (`Dew_lower_bar`) поверх уже
        посчитанного `self.results` — см. докстринг модуля за обоснованием.
        Вызывается только из `calculate()`, один раз, после основного march'а.

        1. Находит первый индекс, где основная ветка уже не bubble (первый
           NaN в `Bubble_bar`) — раньше этой точки ретроградной области
           физически быть не может (см. докстринг модуля). Если такого
           индекса нет (весь диапазон — bubble), делать нечего.
        2. Пробует затравить нижнюю точку через `SaturationPointSSM.
           find_dew_point` на нескольких соседних T подряд начиная с этого
           индекса (первые попытки могут провалиться — ретроградная область
           там ещё слишком узкая) — до `lower_dew_bootstrap_attempts` раз.
        3. От найденной затравки продолжает чистым Ньютоном
           (`SaturationPointNewton.find_dew_point`) до конца диапазона, с той
           же верификацией/бюджетом времени на шаг, что и основная ветка.
        """
        bubble = self.results['Bubble_bar']
        start_idx = next((i for i, v in enumerate(bubble) if pd.isna(v)), None)
        if start_idx is None:
            return

        seed_idx, seed_p = None, None
        attempts = min(self.lower_dew_bootstrap_attempts, len(self.temps_c) - start_idx)
        for j in range(start_idx, start_idx + attempts):
            t_c = self.temps_c[j]
            t_k = float(t_c) + 273.15
            ssm_finder = SaturationPointSSM(self.composition, t_k, self.p_max_bar, self.p_min_bar)
            pt = ssm_finder.find_dew_point(self.p_max_bar)
            if pt is not None and pt.point_type == 'dew':
                seed_idx, seed_p = j, pt.pressure
                logger.info(
                    "Newton: нижняя ретроградная ветка затравлена через SSM на "
                    "T=%.2f °C, P=%.4f бар", t_c, seed_p,
                )
                break

        if seed_idx is None:
            logger.info(
                "Newton: не удалось затравить нижнюю ретроградную ветку через SSM "
                "в первых %d точках после T=%.2f °C", attempts, self.temps_c[start_idx],
            )
            return

        prev_p = seed_p
        slow_streak = 0
        for j in range(seed_idx, len(self.temps_c)):
            t_c = self.temps_c[j]
            t_k = float(t_c) + 273.15
            t0 = time.time()

            finder = SaturationPointNewton(
                self.composition, t_k, self.p_max_bar, self.p_min_bar,
                p_guess_dew=prev_p, tol=self.tol, max_iter=self.max_iter,
            )
            pt = finder.find_dew_point(self.p_max_bar)

            elapsed = time.time() - t0
            slow_streak = slow_streak + 1 if elapsed > self.slow_point_budget_s else 0

            if pt is None:
                logger.info(
                    "Newton: нижняя ретроградная ветка остановлена на T=%.2f °C "
                    "(точка не верифицирована)", t_c,
                )
                break

            upper_here = self.results['Dew_upper_bar'][j]
            if pd.isna(upper_here):
                upper_here = self.results['Bubble_bar'][j]
            if not pd.isna(upper_here) and abs(pt.pressure - upper_here) / upper_here < self.dew_bubble_min_gap_rel:
                logger.info(
                    "Newton: нижняя ретроградная ветка остановлена на T=%.2f °C "
                    "(слилась с верхней — вероятно, конец ретроградной области)", t_c,
                )
                break

            self.results['Dew_lower_bar'][j] = pt.pressure
            prev_p = pt.pressure

            if slow_streak >= self.max_consecutive_slow_points:
                logger.info(
                    "Newton: нижняя ретроградная ветка остановлена на T=%.2f °C "
                    "(%d шагов подряд дольше %.1f сек)",
                    t_c, slow_streak, self.slow_point_budget_s,
                )
                break

    def calculate_parallel(self, n_jobs: int = -1, backend: str = 'loky') -> pd.DataFrame:
        """
        То же самое, что `calculate()`, но диапазон T разбивается на
        непрерывные куски (по числу воркеров), каждый считается независимо со
        своим холодным стартом на первой точке куска — см.
        `_calculate_envelope_chunk_worker`. Если первая точка куска уже
        лежит за критической/ретроградной областью, холодный старт может не
        найтись, и весь кусок вернётся пустым (см. докстринг модуля) — для
        гарантированного покрытия всего диапазона используйте `calculate()`.
        """
        n_workers = max(1, min(effective_n_jobs(n_jobs), len(self.temps_c)))
        chunks = [chunk for chunk in np.array_split(self.temps_c, n_workers) if len(chunk) > 0]

        logger.info(
            "Расчёт фазовой огибающей (Newton, параллельно, %d кусков на %d воркеров): "
            "T=[%.1f, %.1f] °C, шаг %.1f °C",
            len(chunks), n_workers, self.t_min_c, self.t_max_c, self.t_step_c,
        )

        chunk_results = Parallel(n_jobs=n_jobs, backend=backend)(
            delayed(_calculate_envelope_chunk_worker)(
                self.composition, chunk.tolist(), self.p_max_bar, self.p_min_bar,
                self.tol, self.max_iter, self.slow_point_budget_s, self.max_consecutive_slow_points,
            )
            for chunk in chunks
        )

        self.results = {'Temp_C': [], 'Bubble_bar': [], 'Dew_upper_bar': [], 'Dew_lower_bar': []}
        for chunk_result in chunk_results:
            for t_c, primary, secondary in chunk_result:
                secondary = self._filter_near_duplicate_dew(t_c, primary, secondary)
                self._append_result(t_c, primary, secondary)

        return pd.DataFrame(self.results)

    def get_data(self) -> pd.DataFrame:
        return pd.DataFrame(self.results)

    def plot(self, show: bool = True, save_path: str | None = None):
        """P-T диаграмма — тот же вид, что у `PhaseEnvelopeSSM.plot`, три
        физически различные кривые (см. докстринг класса)."""
        df = pd.DataFrame(self.results)
        if df.empty:
            logger.warning("Нет данных для построения — вызовите calculate() перед plot().")
            return None

        fig, ax = plt.subplots(figsize=(10, 6))

        bubble = df.dropna(subset=['Bubble_bar'])
        dew_upper = df.dropna(subset=['Dew_upper_bar'])
        dew_lower = df.dropna(subset=['Dew_lower_bar'])

        if not bubble.empty:
            ax.plot(bubble['Temp_C'], bubble['Bubble_bar'], 'b-o',
                    label='Bubble Point (насыщение)', markersize=4, linewidth=2)
        if not dew_upper.empty:
            ax.plot(dew_upper['Temp_C'], dew_upper['Dew_upper_bar'], 'r-s',
                    label='Dew Point, upper (начало конденсации)', markersize=4, linewidth=2)
        if not dew_lower.empty:
            ax.plot(dew_lower['Temp_C'], dew_lower['Dew_lower_bar'], 'r--^',
                    label='Dew Point, lower (ретроградная)', markersize=4, linewidth=2)

        ax.set_xlabel('Температура, °C')
        ax.set_ylabel('Давление, бар')
        ax.set_title('Фазовая огибающая (метод Ньютона по производным летучести)')
        ax.legend()
        ax.grid(True, linestyle='--', alpha=0.6)
        fig.tight_layout()

        if save_path:
            fig.savefig(save_path, dpi=300, bbox_inches='tight')
            logger.info("График сохранён: %s", save_path)
        if show:
            plt.show()

        return fig
