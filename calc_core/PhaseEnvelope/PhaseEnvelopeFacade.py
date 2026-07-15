"""
Тонкий фасад над каноническими калькуляторами фазовой огибающей.

Привязан к `CompositionalModel` (доступен как `model.phase_envelope`) —
берёт состав из `model.composition`, вызывает нужный calc_core-калькулятор
(`BubblePointCalculator`, `DewPointCalculator`, `CriticalPointCalculator`,
`PhaseEnvelopeGrid`, `PhaseEnvelopeSSM`, `PhaseEnvelopeNewton`) и логирует
каждый вызов в `model.result_store_object` — по тому же паттерну, что уже
использует `CompositionalModel.flash()`.

Единообразие, которое даёт этот слой: у calc_core-калькуляторов разные имена
основного метода (`calculate()`/`calculate_parallel()`/`run_parallel()`) и
разный набор обязательных параметров — здесь это спрятано за одинаковыми по
форме методами `bubble_point`/`dew_point`/`critical_point`/`grid`/`ssm`/`newton`.

Для `bubble_point`/`dew_point`/`critical_point` (одноточечные решения методом
Ньютона) явный сигнал несходимости (`.converged` / `None`) переводится в
`ConvergenceError`, а не молча возвращается наружу. Для `ssm`/`newton` — не
переводится: результат там — таблица по сетке температур, где NaN на части
строк (выше крикондентермы, например) — нормальный физический результат, а не
ошибка расчёта. Для `grid()` — сигнала несходимости нет вовсе (это скан
стабильности, не итеративный солвер).

Валидация входных P/T и общая иерархия исключений (`PVTCalcError` и т.п.) —
не здесь, это отдельная задача (см. `docs/BACKLOG.md`).
"""

import logging

from calc_core.PhaseEnvelope.BubblePointPressure import BubblePointCalculator
from calc_core.PhaseEnvelope.DewPressure import DewPointCalculator
from calc_core.PhaseEnvelope.CriticalPoint import CriticalPointCalculator
from calc_core.PhaseEnvelope.PhaseEnvelopeGrid import PhaseEnvelopeGrid
from calc_core.PhaseEnvelope.PhaseEnvelopeSuccessiveSubstitution import PhaseEnvelopeSSM
from calc_core.PhaseEnvelope.PhaseEnvelopeNewton import PhaseEnvelopeNewton
from calc_core.Utils.Errors import ConvergenceError

logger = logging.getLogger(__name__)


class PhaseEnvelopeFacade:
    """Доступ к каноническим калькуляторам фазовой огибающей через `model.phase_envelope`."""

    def __init__(self, model: 'CompositionalModel'):
        """
        Parameters
        ----------
        model : CompositionalModel
            Родительская модель — состав берётся как `model.composition`,
            результаты логируются в `model.result_store_object`.
        """
        self._model = model

    def bubble_point(self, T: float, **kwargs) -> float:
        """
        Давление точки кипения (Ньютон по производным летучести).

        Parameters
        ----------
        T : float
            Температура, K.
        **kwargs
            Прочие параметры `BubblePointCalculator` (`P_guess`, `K_guess`,
            `max_iter`, `tol`).

        Returns
        -------
        float
            Давление точки кипения, бар.

        Raises
        ------
        ConvergenceError
            Если `BubblePointCalculator` не сошёлся.
        """
        calc = BubblePointCalculator(self._model.composition, T, **kwargs)
        p = calc.calculate()

        self._model.result_store_object.add(
            module='PhaseEnvelope.bubble_point',
            params={'T': T, **kwargs},
            data={'P': p, 'converged': calc.converged},
        )

        if not calc.converged:
            raise ConvergenceError(f'BubblePointCalculator не сошёлся при T={T}')
        return p

    def dew_point(self, T: float, **kwargs) -> float:
        """
        Давление начала конденсации (Ньютон по производным летучести).

        Parameters
        ----------
        T : float
            Температура, K.
        **kwargs
            Прочие параметры `DewPointCalculator` (`dew_point_type`,
            `P_guess`, `K_guess`, `max_iter`, `tol`).

        Returns
        -------
        float
            Давление начала конденсации, бар.

        Raises
        ------
        ConvergenceError
            Если `DewPointCalculator` не сошёлся или вернул `None`.
        """
        calc = DewPointCalculator(self._model.composition, T, **kwargs)
        p = calc.calculate()

        self._model.result_store_object.add(
            module='PhaseEnvelope.dew_point',
            params={'T': T, **kwargs},
            data={'P': p, 'converged': calc.converged},
        )

        if not calc.converged or p is None:
            raise ConvergenceError(f'DewPointCalculator не сошёлся при T={T}')
        return p

    def critical_point(self, **kwargs) -> dict:
        """
        Критическая точка (минимизация зазора bubble/dew).

        Parameters
        ----------
        **kwargs
            Прочие параметры `CriticalPointCalculator` (`verbose`).

        Returns
        -------
        dict
            `{'T_crit', 'P_crit', 'S_v', 'S_l', 'K_v', 'K_l', 'metric',
            'is_trivial_converged'}` — см. `CriticalPointCalculator.calculate()`.

        Raises
        ------
        ConvergenceError
            Если критическая точка не найдена (`calculate()` вернул `None`).
        """
        calc = CriticalPointCalculator(self._model.composition, **kwargs)
        result = calc.calculate()

        self._model.result_store_object.add(
            module='PhaseEnvelope.critical_point',
            params=kwargs,
            data=result,
        )

        if result is None:
            raise ConvergenceError('CriticalPointCalculator не нашёл критическую точку')
        return result

    def grid(self, **kwargs) -> PhaseEnvelopeGrid:
        """
        Сеточный скан стабильности P×T (затравка нижней ветки SSM).

        Parameters
        ----------
        **kwargs
            Параметры `PhaseEnvelopeGrid` (`max_pressure`, `max_temperature`,
            `pressure_points`, `temperature_points`) — все с дефолтами.

        Returns
        -------
        PhaseEnvelopeGrid
            Сам calc-объект после `run_parallel()` — сохраняет `.plot()` и
            прочие атрибуты. В историю (`result_store_object`) кладётся не
            он сам, а снэпшот его массивов на момент вызова.
        """
        calc = PhaseEnvelopeGrid(self._model.composition, **kwargs)
        calc.run_parallel()

        self._model.result_store_object.add(
            module='PhaseEnvelope.grid',
            params=kwargs,
            data={
                'P': list(calc.result_pressure_arr),
                'T': list(calc.result_temperature_arr),
                'stable': list(calc.result_stability_flag_arr),
            },
        )
        return calc

    def ssm(self, t_min_c: float, t_max_c: float, t_step_c: float, p_max_bar: float,
            parallel: bool = True, **kwargs):
        """
        Фазовая огибающая методом последовательных приближений (основная реализация).

        Parameters
        ----------
        t_min_c, t_max_c, t_step_c : float
            Диапазон и шаг температур, °C.
        p_max_bar : float
            Верхняя граница давления, бар.
        parallel : bool, optional
            `True` (по умолчанию) — `calculate_parallel()` (joblib),
            `False` — последовательный `calculate()`.
        **kwargs
            Прочие параметры `PhaseEnvelopeSSM` (`p_min_bar`, `dew_bubble_min_gap_rel`).

        Returns
        -------
        pandas.DataFrame
            Таблица `Temp_C`/`Bubble_bar`/`Dew_upper_bar`/`Dew_lower_bar`.
        """
        calc = PhaseEnvelopeSSM(self._model.composition, t_min_c, t_max_c, t_step_c, p_max_bar, **kwargs)
        df = calc.calculate_parallel() if parallel else calc.calculate()

        self._model.result_store_object.add(
            module='PhaseEnvelope.ssm',
            params={
                't_min_c': t_min_c, 't_max_c': t_max_c, 't_step_c': t_step_c,
                'p_max_bar': p_max_bar, 'parallel': parallel, **kwargs,
            },
            data=df,
        )
        return df

    def newton(self, t_min_c: float, t_max_c: float, t_step_c: float, p_max_bar: float,
               parallel: bool = True, **kwargs):
        """
        Фазовая огибающая методом Ньютона по производным летучести.

        Точнее в отдельной точке, чем `ssm()`, но менее устойчива вблизи
        критической точки.

        Parameters
        ----------
        t_min_c, t_max_c, t_step_c : float
            Диапазон и шаг температур, °C.
        p_max_bar : float
            Верхняя граница давления, бар.
        parallel : bool, optional
            `True` (по умолчанию) — `calculate_parallel()`, `False` — `calculate()`.
        **kwargs
            Прочие параметры `PhaseEnvelopeNewton` (`p_min_bar`, `tol`,
            `max_iter` и др. — см. сигнатуру класса).

        Returns
        -------
        pandas.DataFrame
        """
        calc = PhaseEnvelopeNewton(self._model.composition, t_min_c, t_max_c, t_step_c, p_max_bar, **kwargs)
        df = calc.calculate_parallel() if parallel else calc.calculate()

        self._model.result_store_object.add(
            module='PhaseEnvelope.newton',
            params={
                't_min_c': t_min_c, 't_max_c': t_max_c, 't_step_c': t_step_c,
                'p_max_bar': p_max_bar, 'parallel': parallel, **kwargs,
            },
            data=df,
        )
        return df
