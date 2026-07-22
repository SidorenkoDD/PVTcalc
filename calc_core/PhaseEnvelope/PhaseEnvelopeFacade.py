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

Для `bubble_point`/`dew_point` недостаточно сигнала `.converged`: сырой Newton
может сойтись к математическому корню внутри или вне двухфазной области.
Поэтому фасад дополнительно подтверждает смену устойчивости по обе стороны
давления; неподтверждённый корень переводится в `ConvergenceError`.
`critical_point` пока проверяет только `None` и остаётся отдельной задачей R1.2.
Для `ssm`/`newton` исключение на частичном результате не поднимается: причина
каждого NaN доступна в `DataFrame.attrs['diagnostics']`. Для `grid()` общего
сигнала несходимости пока нет (это набор независимых тестов стабильности).

Граничная валидация входных P/T выполняется здесь; валидация прямых вызовов
низкоуровневых классов остаётся отдельной задачей (см. `docs/BACKLOG.md`).
"""

import logging
from typing import TYPE_CHECKING

from calc_core.PhaseEnvelope.BubblePointPressure import BubblePointCalculator
from calc_core.PhaseEnvelope.CriticalPoint import CriticalPointCalculator
from calc_core.PhaseEnvelope.DewPressure import DewPointCalculator
from calc_core.PhaseEnvelope.PhaseEnvelopeGrid import PhaseEnvelopeGrid
from calc_core.PhaseEnvelope.PhaseEnvelopeNewton import (
    PhaseEnvelopeNewton,
    verify_saturation_boundary,
)
from calc_core.PhaseEnvelope.PhaseEnvelopeSuccessiveSubstitution import PhaseEnvelopeSSM
from calc_core.Utils.Errors import ConvergenceError, InputValidationError, StopIterationError
from calc_core.Utils.Validation import (
    validate_positive_pressure,
    validate_temperature_celsius,
    validate_temperature_kelvin,
)

if TYPE_CHECKING:
    from calc_core.CompositionalModel.CompositionalModel import CompositionalModel

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

    def _composition_copy(self):
        """
        Независимая копия `model.composition` для передачи в calc_core-калькулятор.

        Почти все calc_core-калькуляторы мутируют переданный им `Composition`
        напрямую (`composition.T = ...` — тяжёлый setter, пересчитывает
        EOS-зависимые свойства/BIP, см. `Composition.T`). Без этой копии
        `model.composition` тихо менялся бы при каждом вызове фасада.
        """
        composition = self._model.composition
        return composition.new_composition(composition.composition, deep_copy=True)

    @staticmethod
    def _validate_envelope_range(t_min_c, t_max_c, t_step_c, p_max_bar):
        """Валидация общего для `ssm`/`newton` диапазона огибающей (T в °C, P в бар)."""
        validate_temperature_celsius(t_min_c, name='t_min_c')
        validate_temperature_celsius(t_max_c, name='t_max_c')
        validate_positive_pressure(t_step_c, name='t_step_c')
        validate_positive_pressure(p_max_bar, name='p_max_bar')
        if t_max_c <= t_min_c:
            raise InputValidationError(
                f'Диапазон температур пуст: t_max_c ({t_max_c}) должно быть больше t_min_c ({t_min_c}).')

    def bubble_point(self, T_k: float, **kwargs) -> float:
        """
        Давление точки кипения (Ньютон по производным летучести).

        Parameters
        ----------
        T_k : float
            Температура, **K** (в отличие от `ssm`/`newton`/`grid`, которые
            принимают °C — суффикс в имени параметра именно поэтому).
        **kwargs
            Прочие параметры `BubblePointCalculator` (`P_guess`, `K_guess`,
            `max_iter`, `tol`).

        Returns
        -------
        float
            Давление точки кипения, бар.

        Raises
        ------
        InputValidationError
            Если `T_k` <= 0 K.
        ConvergenceError
            Если калькулятор не сошёлся или корень не подтверждён независимым
            двусторонним тестом стабильности.
        """
        validate_temperature_kelvin(T_k)
        calc = BubblePointCalculator(self._composition_copy(), T_k, **kwargs)
        p = calc.calculate()

        verified_type = None
        verification_error = None
        if calc.converged and p is not None:
            try:
                verified_type = verify_saturation_boundary(
                    calc.composition, T_k, p, nudge_sign=-1.0,
                )
            except StopIterationError as exc:
                verification_error = str(exc)

        self._model.result_store_object.add(
            module='PhaseEnvelope.bubble_point',
            params={'T_k': T_k, **kwargs},
            data={
                'P': p,
                'converged': calc.converged,
                'verified_type': verified_type,
                'verification_error': verification_error,
            },
        )

        if not calc.converged:
            raise ConvergenceError(f'BubblePointCalculator не сошёлся при T={T_k} K')
        if verified_type != 'bubble':
            reason = 'тест стабильности не сошёлся' if verification_error else (
                f'граница не подтверждена как bubble (тип={verified_type})'
            )
            raise ConvergenceError(
                f'BubblePointCalculator вернул непроверенный корень при T={T_k} K: {reason}'
            )
        return p

    def dew_point(self, T_k: float, **kwargs) -> float:
        """
        Давление начала конденсации (Ньютон по производным летучести).

        Parameters
        ----------
        T_k : float
            Температура, **K** (см. замечание про единицы в `bubble_point`).
        **kwargs
            Прочие параметры `DewPointCalculator` (`dew_point_type`,
            `P_guess`, `K_guess`, `max_iter`, `tol`).

        Returns
        -------
        float
            Давление начала конденсации, бар.

        Raises
        ------
        InputValidationError
            Если `T_k` <= 0 K.
        ConvergenceError
            Если калькулятор не сошёлся, вернул `None` или корень не
            подтверждён независимым двусторонним тестом стабильности.
        """
        validate_temperature_kelvin(T_k)
        calc = DewPointCalculator(self._composition_copy(), T_k, **kwargs)
        p = calc.calculate()

        verified_type = None
        verification_error = None
        if calc.converged and p is not None:
            nudge_sign = 1.0 if calc.dew_point_type == 'lower' else -1.0
            try:
                verified_type = verify_saturation_boundary(
                    calc.composition, T_k, p, nudge_sign=nudge_sign,
                )
            except StopIterationError as exc:
                verification_error = str(exc)

        self._model.result_store_object.add(
            module='PhaseEnvelope.dew_point',
            params={'T_k': T_k, **kwargs},
            data={
                'P': p,
                'converged': calc.converged,
                'verified_type': verified_type,
                'verification_error': verification_error,
            },
        )

        if not calc.converged or p is None:
            raise ConvergenceError(f'DewPointCalculator не сошёлся при T={T_k} K')
        if verified_type != 'dew':
            reason = 'тест стабильности не сошёлся' if verification_error else (
                f'граница не подтверждена как dew (тип={verified_type})'
            )
            raise ConvergenceError(
                f'DewPointCalculator вернул непроверенный корень при T={T_k} K: {reason}'
            )
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
            `{'T', 'P', 'T_C', 'P_bubble', 'P_dew', 'gap'}`:
            `T` — температура критической точки, K; `T_C` — она же в °C;
            `P` — давление, бар (полусумма `P_bubble`/`P_dew`);
            `gap` — |P_bubble - P_dew|, **метрика качества сходимости**.
            См. `CriticalPointCalculator.calculate()`.

        Raises
        ------
        ConvergenceError
            Если критическая точка не найдена (`calculate()` вернул `None`).

        Warnings
        --------
        **Результат сейчас недостоверен для многокомпонентных составов.**
        На `KRSNL_PVTSIM` замерен `gap` ≈ 772 бар при координате ~229 °C /
        608 бар — то есть решение разошлось, а `calculate()` всё равно вернул
        словарь, а не `None`. Поэтому GUI намеренно не считает и не рисует
        критическую точку (см. `gui/services/phase_envelope_service.py`).
        **Всегда проверяйте `gap`** перед использованием результата: он должен
        быть малым по сравнению с `P`. Починка солвера — отдельная задача,
        см. `docs/BACKLOG.md`.
        """
        calc = CriticalPointCalculator(self._composition_copy(), **kwargs)
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

        Raises
        ------
        InputValidationError
            Если переданы `max_pressure` <= 0.
        """
        if 'max_pressure' in kwargs:
            validate_positive_pressure(kwargs['max_pressure'], name='max_pressure')
        calc = PhaseEnvelopeGrid(self._composition_copy(), **kwargs)
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

        Raises
        ------
        InputValidationError
            Если `p_max_bar`/`t_step_c` <= 0, температуры ниже абсолютного
            нуля или `t_max_c` <= `t_min_c`.
        """
        self._validate_envelope_range(t_min_c, t_max_c, t_step_c, p_max_bar)
        calc = PhaseEnvelopeSSM(self._composition_copy(), t_min_c, t_max_c, t_step_c, p_max_bar, **kwargs)
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

        Raises
        ------
        InputValidationError
            Если `p_max_bar`/`t_step_c` <= 0, температуры ниже абсолютного
            нуля или `t_max_c` <= `t_min_c`.
        """
        self._validate_envelope_range(t_min_c, t_max_c, t_step_c, p_max_bar)
        calc = PhaseEnvelopeNewton(self._composition_copy(), t_min_c, t_max_c, t_step_c, p_max_bar, **kwargs)
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
