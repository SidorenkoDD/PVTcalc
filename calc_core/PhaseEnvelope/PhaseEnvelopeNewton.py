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

Изначально это списывалось на баг знака в аналитической производной
`dF/dP` (`calc_d_log_phi_i_dp`/"Sign correction" — см. побочную находку в
конце докстринга) — производная действительно была неверна (сверка с
конечными разностями, исправлена в `BrusilovskiyEOS`), но её исправление
СКОРОСТЬ подняло резко (3-7 итераций вместо 20-70 в "спокойной" зоне), а
именно этот эффект (уход к корню внутри области) не убрало. Настоящая
причина — ДРУГАЯ и глубже: `BubblePointCalculator`/`DewPointCalculator` при
КАЖДОМ вызове пересобирают стартовые K-факторы заново по грубой корреляции
Вильсона от переданного P, а не принимают continuation по K с предыдущего
шага — вблизи критической точки состава Вильсон систематически даёт K,
далёкий от истинного равновесного, и Ньютон (даже с верной производной и
точным P) уходит к другому математическому корню уравнения. Добавлены
`K_guess`/`result_K` в оба калькулятора (`BubblePointPressure.py`/
`DewPressure.py`) — march теперь продолжает не только P, но и K между
соседними шагами по T (см. `SaturationPointNewton.last_k` и его
использование в `PhaseEnvelopeNewton.calculate()`), это и есть основной
инструмент, которым ниже "проезжается" трудная зона у критической точки.
Ни то, ни другое не отменяет необходимость двусторонней верификации —
доверять `converged=True` по-прежнему нельзя.

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
совпадение с SSM до 4-5 значащих цифр). Затравка ищется не бисекцией
(`SaturationPointSSM`, дороже), а через ДЕШЁВЫЙ сеточный скан стабильности —
`PhaseEnvelopeGrid.py` (простой P×T скан `TwoPhaseStabilityTest` без
бисекции, ~4.3 сек на 50×97 точек, joblib-параллельно) — `_build_stability_grid`
находит по каждому T-столбцу ВСЕ смены флага стабильности; второй (нижний)
переход в столбце — надёжный дешёвый признак "здесь есть ретроградная зона",
его середина — сразу `P_guess` для `SaturationPointNewton.find_dew_point`
(который сам верифицирует находку — отдельного шага верификации затравки не
нужно). `calculate()` (только последовательный вариант, не
`calculate_parallel()`) после остановки основного march'а вызывает сетку ОДИН
раз (`_build_stability_grid`) и по ней ищет первый T с подтверждённой
ретроградной зоной, дальше продолжает уже чистым Ньютоном continuation'ом по
P и K до конца запрошенного диапазона — с той же верификацией и тем же
бюджетом времени на шаг, что и у основной ветки. `PhaseEnvelopeGrid`
(как и `PhaseEnvelopeSSM`) не редактируется ни строкой — только
импортируется и вызывается как есть.

Для верхней ветки (`_coast_through_upper_dew`) первая версия делала
СИММЕТРИЧНОЕ нижней ветке: одна затравка через `SaturationPointSSM.
find_bubble_point`, потом Ньютон. Не сработала — SSM (бисекция по флагу
стабильности) успешно "находит" точку почти на любом T в трудной зоне
(T≈360-435°C на KRSNL), но `DewPointCalculator`, стартуя оттуда, всё равно
сходился не к границе (см. выше про Вильсон-K) — затравка была не тем
узким местом. Реальным узким местом было отсутствие continuation по K.

С `K_guess` текущая версия (`_coast_through_upper_dew`) устроена иначе — без
SSM вообще: продолжает P И K с ПОСЛЕДНЕЙ уже верифицированной точки
основного march'а (обычно ещё bubble, вплотную к критике) шаг за шагом по
dew-уравнению, но, в отличие от основного march'а, НЕ останавливается на
первой же неверифицированной точке — "проезжает" через несколько таких
шагов подряд (P/K продолжают передаваться дальше, даже когда точка не
записывается в результат), пока уравнение снова не станет хорошо
обусловленным. Проверено эмпирически на составе KRSNL: от конца основного
march'а (T=375°C, bubble) требуется "проехать" через T≈380-395°C без единой
верифицированной точки, прежде чем T=400-425°C снова начинают сходиться и
верифицироваться (согласие с SSM — 0.1-1%). Останавливается, если много
шагов подряд (`max_unverified_run`) не верифицируются ИЛИ Ньютон вовсе не
сходится — то, что осталось непокрытым (включая саму трудную зону, даже
если её "проехали" без записи точек), — по-прежнему честный NaN, покрывается
`PhaseEnvelopeSSM`.

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

ИСТОРИЯ НАХОДКИ (для контекста — сами баги уже исправлены, не в этом файле):
изначально у `BubblePointCalculator` производная `dF/dP` почти на каждой
итерации получалась "не того" знака и принудительно инвертировалась веткой
`if dF_dP > 0: dF_dP = -abs(dF_dP)` ("Sign correction") — расследование
показало, что дело не в этой строке, а в самой аналитической производной
`calc_d_log_phi_i_dp` (`BrusilovskiyEOS._calc_dlogphi_dp_vector`) — она была
неверна (сверка с производной по конечным разностям разошлась и по знаку, и
по величине). Исправлена на численную (центральная разность), "Sign
correction" убрана как ставший ненужным костыль — с этим `BubblePointCalculator`
сходится за 3-7 итераций вместо 20-70. Отдельно от этого — тот самый уход к
корню ВНУТРИ двухфазной области у критической точки (пример с T=390°C выше)
это НЕ объясняло и НЕ исправляло: причина оказалась в отсутствии continuation
по K (см. выше) — исправлено добавлением `K_guess`/`result_K` в оба
калькулятора.
"""

import logging
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from joblib import Parallel, delayed, effective_n_jobs

from calc_core.Composition.Composition import Composition
from calc_core.PhaseEnvelope.BubblePointPressure import BubblePointCalculator
from calc_core.PhaseEnvelope.DewPressure import DewPointCalculator
from calc_core.PhaseEnvelope.Diagnostics import (
    NOT_PRESENT,
    NOT_SEARCHED,
    OK,
    OUT_OF_RANGE,
    SKIPPED_AFTER_FAILURE,
    SLOW_POINT_LIMIT,
    SOLVER_NONCONVERGENCE,
    STABILITY_NONCONVERGENCE,
    VERIFICATION_FAILED,
    EnvelopeDiagnostic,
    attach_diagnostics,
)
from calc_core.PhaseEnvelope.PhaseEnvelopeGrid import PhaseEnvelopeGrid
from calc_core.PhaseEnvelope.PhaseEnvelopeSuccessiveSubstitution import SaturationPoint
from calc_core.PhaseStability.TwoPhaseStabilityTest import TwoPhaseStabilityTest
from calc_core.Utils.Errors import StopIterationError

logger = logging.getLogger(__name__)


def verify_saturation_boundary(
    composition: Composition,
    t_k: float,
    p_bar: float,
    nudge_sign: float,
    verify_nudge_rel: float = 1e-3,
) -> str | None:
    """Проверяет, что сырой корень действительно разделяет две фазы и одну.

    ``nudge_sign`` указывает сторону двухфазной области: ``-1`` для верхней
    границы (bubble/dew-upper), ``+1`` для нижней ретроградной dew-границы.
    Возвращает физический тип границы либо ``None``, если обе стороны имеют
    одинаковый флаг стабильности. Несходимость теста стабильности не
    маскируется и передаётся вызывающему коду как ``StopIterationError``.
    """
    p_inside = p_bar * (1.0 + nudge_sign * 10.0 * verify_nudge_rel)
    p_outside = p_bar * (1.0 - nudge_sign * 10.0 * verify_nudge_rel)

    stab_inside = TwoPhaseStabilityTest(composition, p_inside, t_k)
    stab_inside.calculate_phase_stability()
    if stab_inside.stable:
        return None

    stab_outside = TwoPhaseStabilityTest(composition, p_outside, t_k)
    stab_outside.calculate_phase_stability()
    if not stab_outside.stable:
        return None

    return 'bubble' if (stab_inside.S_v - 1.0) >= (stab_inside.S_l - 1.0) else 'dew'


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
        p_guess_dew_upper: float | None = None,
        p_guess_dew_lower: float | None = None,
        k_guess_bubble: dict | None = None,
        k_guess_dew_upper: dict | None = None,
        k_guess_dew_lower: dict | None = None,
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
        # dew_upper — guess для dew-ветки ВНУТРИ `find_bubble_point` (fallback
        # первичного/верхнего поиска); dew_lower — guess для `find_dew_point`
        # (нижняя ретроградная точка). Раньше это был ОДИН и тот же атрибут
        # (`p_guess_dew`) — из-за чего верхний dew-fallback внутри march'а по
        # ошибке получал continuation-guess от НИЖНЕЙ ветки (физически другое
        # давление, ~10x меньше) вместо своего собственного — разделены.
        self.p_guess_dew_upper = p_guess_dew_upper
        self.p_guess_dew_lower = p_guess_dew_lower
        # K_guess — continuation по K-факторам с предыдущего сошедшегося шага
        # по T (не только по P). Без этого `BubblePointCalculator`/
        # `DewPointCalculator` при КАЖДОМ вызове пересобирают K заново по
        # грубой корреляции Вильсона от текущего P — вблизи критической точки
        # состава это может увести Ньютон к другому математическому корню
        # даже от точного P_guess (см. докстринг модуля). `last_k` — K
        # последней УСПЕШНО решённой точки (какого бы типа она ни была),
        # чтобы вызывающий код (march в `PhaseEnvelopeNewton`) мог передать
        # его дальше как continuation для следующего T.
        self.k_guess_bubble = k_guess_bubble
        self.k_guess_dew_upper = k_guess_dew_upper
        self.k_guess_dew_lower = k_guess_dew_lower
        self.last_k: dict | None = None
        self.last_primary_status = NOT_SEARCHED
        self.last_secondary_status = NOT_SEARCHED
        self.last_verification_status = NOT_SEARCHED
        self.branch_hint = branch_hint if branch_hint in ('bubble', 'dew') else 'bubble'
        self.tol = tol
        self.max_iter = max_iter
        self.verify_nudge_rel = verify_nudge_rel

    def _solve_bubble(self, guess: float | None, k_guess: dict | None) -> tuple[float | None, bool, dict | None]:
        calc = BubblePointCalculator(
            self.composition, self.t_k, P_guess=guess, K_guess=k_guess,
            tol=self.tol, max_iter=self.max_iter,
        )
        p = calc.calculate()
        return p, calc.converged, calc.result_K

    def _solve_dew(self, guess: float | None, k_guess: dict | None, dew_point_type: str) -> tuple[float | None, bool, dict | None]:
        calc = DewPointCalculator(
            self.composition, self.t_k, dew_point_type=dew_point_type,
            P_guess=guess, K_guess=k_guess, tol=self.tol, max_iter=self.max_iter,
        )
        p = calc.calculate()
        return p, calc.converged, calc.result_K

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
        try:
            label = verify_saturation_boundary(
                self.composition, self.t_k, p, nudge_sign, self.verify_nudge_rel,
            )
        except StopIterationError:
            self.last_verification_status = STABILITY_NONCONVERGENCE
            return None
        self.last_verification_status = OK if label is not None else VERIFICATION_FAILED
        return label

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
        self.last_primary_status = NOT_SEARCHED
        attempted_statuses: list[str] = []
        order = ['bubble', 'dew'] if self.branch_hint == 'bubble' else ['dew', 'bubble']

        for kind in order:
            if kind == 'bubble':
                p, converged, result_k = self._solve_bubble(self.p_guess_bubble, self.k_guess_bubble)
            else:
                p, converged, result_k = self._solve_dew(self.p_guess_dew_upper, self.k_guess_dew_upper, 'upper')

            if not converged or p is None:
                attempted_statuses.append(SOLVER_NONCONVERGENCE)
                continue
            if not (self.p_min_bar < p < self.p_max_bar):
                attempted_statuses.append(OUT_OF_RANGE)
                continue

            label = self._verify(p, nudge_sign=-1.0)
            if label != kind:
                attempted_statuses.append(self.last_verification_status)
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
            self.last_k = result_k
            self.last_primary_status = OK
            return SaturationPoint(pressure=float(p), point_type=label)

        for status in (STABILITY_NONCONVERGENCE, VERIFICATION_FAILED, OUT_OF_RANGE, SOLVER_NONCONVERGENCE):
            if status in attempted_statuses:
                self.last_primary_status = status
                break
        logger.info("Newton: верхняя точка насыщения не найдена/не верифицирована (T=%.2f K)", self.t_k)
        return None

    def find_dew_point(self, p_upper_bound: float) -> SaturationPoint | None:
        """
        Поиск второй (нижней, ретроградной) точки насыщения — только
        dew-уравнение (`dew_point_type='lower'`), только по continuation-guess
        (`p_guess_dew_lower`). Без явного guess практически не сходится к
        верному корню (см. докстринг модуля) — в отличие от `SaturationPointSSM`,
        тут нет бисекции/скана, которые могли бы обнаружить эту точку "вслепую".
        """
        self.last_secondary_status = NOT_SEARCHED
        if self.p_guess_dew_lower is None:
            logger.debug(
                "Newton: нижняя точка насыщения не ищется без guess (T=%.2f K) — "
                "у метода нет надёжного холодного старта для этой ветки",
                self.t_k,
            )
            return None

        p, converged, result_k = self._solve_dew(self.p_guess_dew_lower, self.k_guess_dew_lower, 'lower')
        if not converged or p is None:
            self.last_secondary_status = SOLVER_NONCONVERGENCE
            return None
        if not (self.p_min_bar < p < min(p_upper_bound, self.p_max_bar)):
            self.last_secondary_status = OUT_OF_RANGE
            return None

        label = self._verify(p, nudge_sign=1.0)
        if label != 'dew':
            self.last_secondary_status = self.last_verification_status
            logger.debug(
                "Newton: нижняя точка P=%.4f бар (T=%.2f K) не прошла "
                "верификацию (метка=%s) — отбраковка", p, self.t_k, label,
            )
            return None

        logger.info("Newton: нижняя точка насыщения P=%.4f бар (T=%.2f K)", p, self.t_k)
        self.last_k = result_k
        self.last_secondary_status = OK
        return SaturationPoint(pressure=float(p), point_type='dew')


def _calculate_envelope_chunk_worker(
    composition: Composition, t_c_chunk: list[float], p_max_bar: float, p_min_bar: float,
    tol: float, max_iter: int, slow_point_budget_s: float, max_consecutive_slow_points: int,
) -> list[tuple[float, SaturationPoint | None, SaturationPoint | None, str, str]]:
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
    prev_p_dew_upper: float | None = None
    prev_p_dew_lower: float | None = None
    prev_kb: dict | None = None
    prev_k_dew_upper: dict | None = None
    prev_k_dew_lower: dict | None = None
    prev_type = 'bubble'
    slow_streak = 0

    for i, t_c in enumerate(t_c_chunk):
        t_k = float(t_c) + 273.15
        t0 = time.time()

        finder = SaturationPointNewton(
            comp, t_k, p_max_bar, p_min_bar,
            p_guess_bubble=prev_pb, p_guess_dew_upper=prev_p_dew_upper,
            p_guess_dew_lower=prev_p_dew_lower,
            k_guess_bubble=prev_kb, k_guess_dew_upper=prev_k_dew_upper,
            k_guess_dew_lower=prev_k_dew_lower,
            branch_hint=prev_type, tol=tol, max_iter=max_iter,
        )
        primary = finder.find_bubble_point()
        primary_k = finder.last_k
        p_upper_for_dew = primary.pressure if primary is not None else p_max_bar
        secondary = finder.find_dew_point(p_upper_for_dew)
        secondary_k = finder.last_k if secondary is not None else None

        elapsed = time.time() - t0
        slow_streak = slow_streak + 1 if elapsed > slow_point_budget_s else 0

        if primary is None or slow_streak >= max_consecutive_slow_points:
            primary_status = (
                finder.last_primary_status if primary is None else SLOW_POINT_LIMIT
            )
            results.append((
                float(t_c), None, None,
                primary_status, finder.last_secondary_status,
            ))
            results.extend(
                (float(rest_t), None, None, SKIPPED_AFTER_FAILURE, SKIPPED_AFTER_FAILURE)
                for rest_t in t_c_chunk[i + 1:]
            )
            break

        results.append((
            float(t_c), primary, secondary,
            finder.last_primary_status, finder.last_secondary_status,
        ))
        prev_type = primary.point_type
        if primary.point_type == 'bubble':
            prev_pb, prev_kb = primary.pressure, primary_k
        else:
            prev_p_dew_upper, prev_k_dew_upper = primary.pressure, primary_k
        if secondary is not None:
            prev_p_dew_lower, prev_k_dew_lower = secondary.pressure, secondary_k

    return results


class PhaseEnvelopeNewton:
    """
    Фазовая огибающая методом Ньютона по производным летучести — второй,
    независимый от `PhaseEnvelopeSSM` способ построить огибающую (см.
    докстринг модуля за мотивацией/ограничениями). Тот же интерфейс и формат
    результата, что у `PhaseEnvelopeSSM` (для прямого сравнения/наложения на
    графике). Основной march (bubble/dew upper с continuation по P и K)
    останавливается у критической точки состава — `calculate()` (не
    `calculate_parallel()`) следом достраивает обе оставшиеся ветки отдельными
    проходами: нижнюю ретроградную (`Dew_lower_bar`) — через одну затравку из
    дешёвого сеточного скана стабильности (`_build_stability_grid`,
    `PhaseEnvelopeGrid`) и дальше continuation Ньютоном по P и K (см.
    `_bootstrap_and_march_lower_dew`); верхнюю (`Dew_upper_bar`) сразу за
    критикой — "проездом" через несколько неверифицированных шагов подряд с
    последней уже верифицированной точки основного march'а (см.
    `_coast_through_upper_dew`). Ни один из этих проходов не использует
    `SaturationPointSSM`/бисекцию. Обе достройки покрывают то, что могут, но
    не гарантируют полного покрытия до самого конца запрошенного диапазона —
    для гарантии в критической/ретроградной зоне остаётся `PhaseEnvelopeSSM`.
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
        grid_pressure_points: int = 50,
        grid_temperature_points: int | None = None,
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
        # Затравка нижней ретроградной ветки через сеточный скан стабильности
        # (см. докстринг модуля и `_build_stability_grid`) — запускается только
        # в calculate() (не в calculate_parallel()), только один раз за march,
        # и только если основная ветка вообще покинула bubble-область в
        # запрошенном диапазоне T.
        self.bootstrap_lower_dew = bootstrap_lower_dew
        # Разрешение сетки `PhaseEnvelopeGrid` для затравки — по P
        # фиксированное, по T по умолчанию совпадает с уже запрошенным шагом
        # march'а (`len(self.temps_c)`), можно переопределить явно (например,
        # для другого состава/диапазона может понадобиться другое разрешение).
        self.grid_pressure_points = grid_pressure_points
        self.grid_temperature_points = grid_temperature_points

        self.temps_c = np.arange(t_min_c, t_max_c + t_step_c / 2.0, t_step_c)
        self.results = {'Temp_C': [], 'Bubble_bar': [], 'Dew_upper_bar': [], 'Dew_lower_bar': []}
        self.diagnostics: list[EnvelopeDiagnostic] = []

    def _frame(self) -> pd.DataFrame:
        return attach_diagnostics(
            pd.DataFrame(self.results), self.diagnostics,
            method="newton", fallback="ssm",
        )

    def _append_diagnostic(self, t_c: float, primary_status: str, secondary_status: str) -> None:
        self.diagnostics.append({
            'Temp_C': float(t_c),
            'primary_status': primary_status,
            'secondary_status': secondary_status,
        })

    def _set_diagnostic_status(
        self,
        index: int,
        *,
        primary_status: str | None = None,
        secondary_status: str | None = None,
    ) -> None:
        if primary_status is not None:
            self.diagnostics[index]['primary_status'] = primary_status
        if secondary_status is not None:
            self.diagnostics[index]['secondary_status'] = secondary_status

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

    def _append_nan_rows(self, temps_c, status: str = SKIPPED_AFTER_FAILURE) -> None:
        for t_c in temps_c:
            self.results['Temp_C'].append(float(t_c))
            self.results['Bubble_bar'].append(np.nan)
            self.results['Dew_upper_bar'].append(np.nan)
            self.results['Dew_lower_bar'].append(np.nan)
            self._append_diagnostic(t_c, status, status)

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
        self.diagnostics = []

        prev_pb: float | None = None
        prev_p_dew_upper: float | None = None
        prev_p_dew_lower: float | None = None
        prev_kb: dict | None = None
        prev_k_dew_upper: dict | None = None
        prev_k_dew_lower: dict | None = None
        prev_type = 'bubble'
        slow_streak = 0
        last_primary_p: float | None = None
        last_primary_k: dict | None = None

        for i, t_c in enumerate(self.temps_c):
            t_k = float(t_c) + 273.15
            t0 = time.time()

            finder = SaturationPointNewton(
                self.composition, t_k, self.p_max_bar, self.p_min_bar,
                p_guess_bubble=prev_pb, p_guess_dew_upper=prev_p_dew_upper,
                p_guess_dew_lower=prev_p_dew_lower,
                k_guess_bubble=prev_kb, k_guess_dew_upper=prev_k_dew_upper,
                k_guess_dew_lower=prev_k_dew_lower,
                branch_hint=prev_type, tol=self.tol, max_iter=self.max_iter,
            )
            primary = finder.find_bubble_point()
            primary_k = finder.last_k
            p_upper_for_dew = primary.pressure if primary is not None else self.p_max_bar
            secondary = finder.find_dew_point(p_upper_for_dew)
            secondary_k = finder.last_k if secondary is not None else None
            secondary = self._filter_near_duplicate_dew(t_c, primary, secondary)
            secondary_status = finder.last_secondary_status if secondary is not None else (
                NOT_PRESENT
                if finder.last_secondary_status == OK
                else finder.last_secondary_status
            )

            elapsed = time.time() - t0
            slow_streak = slow_streak + 1 if elapsed > self.slow_point_budget_s else 0

            if primary is None:
                logger.info(
                    "Newton march остановлен на T=%.2f °C (точка не верифицирована) — "
                    "вероятно, окрестность критической точки состава; для этой зоны "
                    "используйте PhaseEnvelopeSSM", t_c,
                )
                self._append_result(t_c, None, None)
                self._append_diagnostic(
                    t_c, finder.last_primary_status, secondary_status,
                )
                self._append_nan_rows(self.temps_c[i + 1:])
                break

            if slow_streak >= self.max_consecutive_slow_points:
                logger.info(
                    "Newton march остановлен на T=%.2f °C (%d шагов подряд дольше "
                    "%.1f сек) — вероятно, окрестность критической точки состава; "
                    "для этой зоны используйте PhaseEnvelopeSSM",
                    t_c, slow_streak, self.slow_point_budget_s,
                )
                self._append_result(t_c, primary, secondary)
                self._append_diagnostic(t_c, OK, secondary_status)
                self._append_nan_rows(self.temps_c[i + 1:], SLOW_POINT_LIMIT)
                break

            self._append_result(t_c, primary, secondary)
            self._append_diagnostic(t_c, finder.last_primary_status, secondary_status)
            prev_type = primary.point_type
            if primary.point_type == 'bubble':
                prev_pb, prev_kb = primary.pressure, primary_k
            else:
                prev_p_dew_upper, prev_k_dew_upper = primary.pressure, primary_k
            if secondary is not None:
                prev_p_dew_lower, prev_k_dew_lower = secondary.pressure, secondary_k
            # Последняя верифицированная точка march'а (любого типа) — лучшая
            # доступная затравка для "проезда" через трудную зону вокруг
            # критической точки состава (см. `_coast_through_upper_dew`):
            # точнее, чем затравка через SSM с ослабленным допуском.
            last_primary_p, last_primary_k = primary.pressure, primary_k

        if self.bootstrap_lower_dew:
            # Сначала верхняя ветка (Bubble/Dew_upper) — нижняя ветка сверяется
            # с ней при отбраковке дублей (см. ниже), поэтому порядок важен.
            self._coast_through_upper_dew(last_primary_p, last_primary_k)
            grid = self._build_stability_grid()
            self._bootstrap_and_march_lower_dew(grid)

        return self._frame()

    def _build_stability_grid(self) -> dict[float, list[tuple[float, float]]]:
        """
        Грубая (но дешёвая) карта переходов стабильности по сетке P×T —
        затравка для `_bootstrap_and_march_lower_dew`. Использует уже
        существующий `PhaseEnvelopeGrid` (простой скан стабильности
        по решётке, без бисекции — не путать с `PhaseEnvelopeSSM`, тот файл
        тоже не редактируется, только импортируется и вызывается как есть).

        Численно проверено (состав KRSNL, 50×97 точек, полный диапазон
        10-480°C): скан занимает ~4.3 сек и КОРРЕКТНО обнаруживает структуру
        переходов флага стабильности по каждому T-столбцу — 1 переход ниже
        критической точки состава (обычная bubble-граница), 2 перехода в
        ретроградной зоне (нижняя и верхняя граница) и 0 переходов за
        крикондентермом (однофазный газ при любом P — тоже физически верно).
        Именно наличие ВТОРОГО (нижнего) перехода на T — надёжный, дешёвый
        признак "здесь есть ретроградная область", без которого
        `_bootstrap_and_march_lower_dew` раньше приходилось узнавать только
        через дорогую бисекцию `SaturationPointSSM.find_dew_point` (убрана).

        Разрешение по T `PhaseEnvelopeGrid` всегда строит от 0°C
        (`np.linspace(273.15, max_temperature+273.15, ...)`), не от
        `t_min_c` — поэтому её T-сетка в общем случае не совпадает 1:1 с
        `self.temps_c`; сопоставление делается по ближайшему T при
        использовании результата (см. `_bootstrap_and_march_lower_dew`).

        Returns
        -------
        dict[float, list[tuple[float, float]]]
            `{T_c (из сетки PhaseEnvelopeGrid): [(p_lo, p_hi), ...]}`
            — по одной паре `(p_lo, p_hi)` вокруг каждой найденной смены флага
            стабильности на этом T, отсортированные по возрастанию P (первый
            элемент — самый нижний переход, если их несколько).
        """
        grid = PhaseEnvelopeGrid(
            self.composition,
            max_pressure=self.p_max_bar,
            max_temperature=self.t_max_c,
            pressure_points=self.grid_pressure_points,
            temperature_points=self.grid_temperature_points or len(self.temps_c),
        )
        grid.run_parallel()

        p_arr = np.array(grid.result_pressure_arr)
        t_arr = np.array(grid.result_temperature_arr)
        f_arr = np.array(grid.result_stability_flag_arr)

        transitions_by_t: dict[float, list[tuple[float, float]]] = {}
        for t_c in np.unique(t_arr):
            mask = np.isclose(t_arr, t_c)
            p_col = p_arr[mask]
            f_col = f_arr[mask]
            order = np.argsort(p_col)
            p_sorted = p_col[order]
            f_sorted = f_col[order]

            transitions = [
                (float(p_sorted[i]), float(p_sorted[i + 1]))
                for i in range(len(f_sorted) - 1)
                if f_sorted[i] != f_sorted[i + 1]
            ]
            if transitions:
                transitions_by_t[float(t_c)] = transitions

        return transitions_by_t

    def _coast_through_upper_dew(
        self, seed_p: float | None, seed_k: dict | None, max_unverified_run: int = 12,
    ) -> None:
        """
        Пробует "проехать" верхнюю ветку (`Dew_upper_bar`) через трудную зону
        вокруг критической точки состава, где основной march остановился —
        не отдельными затравками через SSM (как было раньше — см. историю
        находки ниже), а continuation'ом P И K с ПОСЛЕДНЕЙ УЖЕ верифицированной
        точки march'а (`seed_p`/`seed_k`), продолжая уравнение "dew" по T шаг
        за шагом БЕЗ остановки на первой же неверифицированной точке.

        Смена философии по сравнению с первой версией этого метода: раньше
        continuation обрывался сразу, как только точка не проходила
        верификацию. Обнаружено эмпирически (ручной прогон вне модуля): если
        ПРОДОЛЖАТЬ continuation через несколько таких неподтверждённых шагов
        подряд (перенося P и K дальше, даже не записывая их в результат), то
        через какое-то число шагов T уравнение "dew" снова становится хорошо
        обусловленным, и точки снова начинают проходить верификацию — на
        составе KRSNL удалось дойти от T=375°C (bubble, конец основного
        march'а) до T=425°C (dew, верифицировано, совпадает с SSM до 0.2-0.5%)
        за 2-3 итерации Ньютона на шаг. Причина, по которой первая версия
        (одна SSM-затравка сразу за разрывом, без "проезда") не работала:
        сам `DewPointCalculator` пересобирает K заново по Вильсону при
        КАЖДОМ вызове без continuation — SSM это ограничение не имеет
        (бисекция по флагу стабильности, не Ньютон), но найти НАСТОЯЩУЮ
        границу от SSM-затравки `DewPointCalculator` всё равно не мог — не
        потому что затравка плохая, а потому что continuation'а по K не было
        вообще (только по P). Теперь есть (`K_guess`/`result_K` в
        `BubblePointCalculator`/`DewPointCalculator`).

        Только ВЕРИФИЦИРОВАННЫЕ точки записываются в `Dew_upper_bar` —
        непройденные остаются NaN (честно), но их P/K всё равно используются
        как continuation-основа для следующего шага (в этом и смысл "проезда").
        Обрывается, если `max_unverified_run` шагов подряд не сходятся вообще
        (Newton не сошёлся) — за пределами этого счётчика зона считается
        непроходимой в разумное время, и `PhaseEnvelopeSSM` — единственный
        вариант для неё.
        """
        if seed_p is None or seed_k is None:
            return

        dew_upper = self.results['Dew_upper_bar']
        bubble = self.results['Bubble_bar']
        start_idx = next(
            (i for i, (b, d) in enumerate(zip(bubble, dew_upper)) if pd.isna(b) and pd.isna(d)),
            None,
        )
        if start_idx is None:
            return  # весь диапазон уже покрыт

        prev_p, prev_k = seed_p, seed_k
        unverified_run = 0

        for j in range(start_idx, len(self.temps_c)):
            t_c = self.temps_c[j]
            t_k = float(t_c) + 273.15

            calc = DewPointCalculator(
                self.composition, t_k, dew_point_type='upper',
                P_guess=prev_p, K_guess=prev_k, tol=self.tol, max_iter=self.max_iter,
            )
            p = calc.calculate()

            if p is None or not calc.converged or not (self.p_min_bar < p < self.p_max_bar):
                failure_status = (
                    SOLVER_NONCONVERGENCE
                    if p is None or not calc.converged
                    else OUT_OF_RANGE
                )
                self._set_diagnostic_status(j, primary_status=failure_status)
                unverified_run += 1
                logger.debug(
                    "Newton (coast-through): не сошлось на T=%.2f °C (%d/%d "
                    "неудач подряд)", t_c, unverified_run, max_unverified_run,
                )
                if unverified_run > max_unverified_run:
                    logger.info(
                        "Newton: верхняя ветка (coast-through) остановлена на "
                        "T=%.2f °C — %d шагов подряд не сходятся вообще",
                        t_c, unverified_run,
                    )
                    return
                continue  # продолжаем со СТАРОЙ (prev_p, prev_k) основы

            verify_finder = SaturationPointNewton(self.composition, t_k, self.p_max_bar, self.p_min_bar)
            label = verify_finder._verify(p, nudge_sign=-1.0)

            # Продолжаем ВСЕГДА (даже без верификации) — в этом и есть смысл
            # "проезда": сам факт сходимости Ньютона (даже к нефизичному
            # корню) обычно ближе к истинной равновесной траектории K(T), чем
            # застревание на старой точке.
            prev_p, prev_k = p, calc.result_K

            if label == 'dew':
                dew_upper[j] = p
                self._set_diagnostic_status(j, primary_status=OK)
                unverified_run = 0
                logger.info(
                    "Newton: верхняя ветка (coast-through) — верифицированная "
                    "точка на T=%.2f °C, P=%.4f бар", t_c, p,
                )
            else:
                self._set_diagnostic_status(
                    j, primary_status=verify_finder.last_verification_status,
                )
                unverified_run += 1
                if unverified_run > max_unverified_run:
                    logger.info(
                        "Newton: верхняя ветка (coast-through) остановлена на "
                        "T=%.2f °C — %d шагов подряд не верифицируются",
                        t_c, unverified_run,
                    )
                    return

    def _bootstrap_and_march_lower_dew(self, grid: dict[float, list[tuple[float, float]]]) -> None:
        """
        Достраивает нижнюю ретроградную dew-ветку (`Dew_lower_bar`) поверх уже
        посчитанного `self.results` — см. докстринг модуля за обоснованием.
        Вызывается только из `calculate()`, один раз, после основного march'а.

        1. Находит первый индекс, где основная ветка уже не bubble (первый
           NaN в `Bubble_bar`) — раньше этой точки ретроградной области
           физически быть не может (см. докстринг модуля). Если такого
           индекса нет (весь диапазон — bubble), делать нечего.
        2. Идёт по T начиная с этого индекса и ищет первый T, для которого
           `grid` (см. `_build_stability_grid`) показывает ВТОРОЙ (нижний)
           переход стабильности — то есть достоверный, дешёвый признак, что
           здесь физически есть ретроградная зона. Середина этого перехода —
           `P_guess` для `SaturationPointNewton.find_dew_point` (холодный
           старт по K, но не по P — уже не Вильсон, а точка из сетки).
        3. От найденной (и уже верифицированной — `find_dew_point` сам это
           делает) затравки продолжает чистым Ньютоном continuation'ом по P
           и K до конца диапазона, с той же верификацией/бюджетом времени на
           шаг, что и основная ветка.
        """
        bubble = self.results['Bubble_bar']
        start_idx = next((i for i, v in enumerate(bubble) if pd.isna(v)), None)
        if start_idx is None:
            return

        if not grid:
            logger.info("Newton: сеточная затравка нижней ветки недоступна — пустая сетка стабильности")
            return

        grid_temps = np.array(sorted(grid.keys()))

        seed_idx, seed_p, seed_k = None, None, None
        for j in range(start_idx, len(self.temps_c)):
            t_c = self.temps_c[j]
            nearest_t = float(grid_temps[np.argmin(np.abs(grid_temps - t_c))])
            transitions = grid[nearest_t]
            if len(transitions) < 2:
                continue  # по сетке здесь ретроградной зоны нет (или мимо разрешения)

            p_lo, p_hi = transitions[0]  # самый нижний переход в этом T-столбце
            guess = (p_lo + p_hi) / 2.0

            t_k = float(t_c) + 273.15
            finder = SaturationPointNewton(
                self.composition, t_k, self.p_max_bar, self.p_min_bar,
                p_guess_dew_lower=guess, tol=self.tol, max_iter=self.max_iter,
            )
            pt = finder.find_dew_point(self.p_max_bar)
            if pt is not None:
                seed_idx, seed_p, seed_k = j, pt.pressure, finder.last_k
                self.results['Dew_lower_bar'][j] = seed_p
                self._set_diagnostic_status(j, secondary_status=OK)
                logger.info(
                    "Newton: нижняя ретроградная ветка затравлена по сетке "
                    "стабильности на T=%.2f °C, P=%.4f бар (сетка: %.2f бар)",
                    t_c, seed_p, guess,
                )
                break
            self._set_diagnostic_status(
                j, secondary_status=finder.last_secondary_status,
            )

        if seed_idx is None:
            logger.info(
                "Newton: сетка не показала подтверждённой ретроградной зоны "
                "в запрошенном диапазоне T (после T=%.2f °C)", self.temps_c[start_idx],
            )
            return

        prev_p, prev_k = seed_p, seed_k
        slow_streak = 0
        for j in range(seed_idx + 1, len(self.temps_c)):
            t_c = self.temps_c[j]
            t_k = float(t_c) + 273.15
            t0 = time.time()

            finder = SaturationPointNewton(
                self.composition, t_k, self.p_max_bar, self.p_min_bar,
                p_guess_dew_lower=prev_p, k_guess_dew_lower=prev_k,
                tol=self.tol, max_iter=self.max_iter,
            )
            pt = finder.find_dew_point(self.p_max_bar)

            elapsed = time.time() - t0
            slow_streak = slow_streak + 1 if elapsed > self.slow_point_budget_s else 0

            if pt is None:
                self._set_diagnostic_status(
                    j, secondary_status=finder.last_secondary_status,
                )
                logger.info(
                    "Newton: нижняя ретроградная ветка остановлена на T=%.2f °C "
                    "(точка не верифицирована)", t_c,
                )
                break

            upper_here = self.results['Dew_upper_bar'][j]
            if pd.isna(upper_here):
                upper_here = self.results['Bubble_bar'][j]
            if not pd.isna(upper_here) and abs(pt.pressure - upper_here) / upper_here < self.dew_bubble_min_gap_rel:
                self._set_diagnostic_status(j, secondary_status=NOT_PRESENT)
                logger.info(
                    "Newton: нижняя ретроградная ветка остановлена на T=%.2f °C "
                    "(слилась с верхней — вероятно, конец ретроградной области)", t_c,
                )
                break

            self.results['Dew_lower_bar'][j] = pt.pressure
            self._set_diagnostic_status(j, secondary_status=OK)
            prev_p, prev_k = pt.pressure, finder.last_k

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
        self.diagnostics = []
        for chunk_result in chunk_results:
            for t_c, primary, secondary, primary_status, secondary_status in chunk_result:
                secondary = self._filter_near_duplicate_dew(t_c, primary, secondary)
                if secondary is None and secondary_status == OK:
                    secondary_status = NOT_PRESENT
                self._append_result(t_c, primary, secondary)
                self._append_diagnostic(t_c, primary_status, secondary_status)

        return self._frame()

    def get_data(self) -> pd.DataFrame:
        return self._frame()

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
