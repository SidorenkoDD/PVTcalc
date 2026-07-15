"""
Валидаторы входных данных для сервисного слоя (фасады `CompositionalModel`).

Чистые функции, бросающие `InputValidationError` при некорректном входе —
единый источник правды, переиспользуемый `CompositionalModel.flash`,
`PhaseEnvelopeFacade` и `ExperimentsFacade`. Валидация вызывается на границе
фасада, до создания calc_core-объектов, чтобы отрицательное/нулевое давление,
нефизичная температура или ненормированный состав давали понятную ошибку, а
не NaN/мусор из numpy/scipy.

Температура намеренно разведена на две функции (`_kelvin`/`_celsius`), т.к.
разные методы фасадов принимают её в разных единицах (см. их сигнатуры) —
вызывающий выбирает подходящую.
"""

from typing import Iterable, Union

from calc_core.Utils.Errors import InputValidationError

# Абсолютный ноль в °C — физическая нижняя граница температуры (эквивалент T_K > 0).
_ABSOLUTE_ZERO_C = -273.15


def validate_positive_pressure(p: Union[float, Iterable[float]], name: str = 'давление') -> None:
    """
    Проверяет, что давление(я) строго положительно.

    Parameters
    ----------
    p : float | Iterable[float]
        Скаляр или последовательность ступеней давления, бар.
    name : str
        Имя параметра для сообщения об ошибке.

    Raises
    ------
    InputValidationError
        Если значение (или любой элемент) <= 0.
    """
    values = p if isinstance(p, Iterable) else [p]
    for value in values:
        if value <= 0.0:
            raise InputValidationError(
                f'{name} должно быть больше 0 (бар). Передано: {value}')


def validate_temperature_kelvin(t: float, name: str = 'температура') -> None:
    """
    Проверяет, что температура (в Кельвинах) строго положительна.

    Raises
    ------
    InputValidationError
        Если `t` <= 0.
    """
    if t <= 0.0:
        raise InputValidationError(
            f'{name} должна быть больше 0 K. Передано: {t}')


def validate_temperature_celsius(t: float, name: str = 'температура') -> None:
    """
    Проверяет, что температура (в °C) выше абсолютного нуля (эквивалент T_K > 0).

    Raises
    ------
    InputValidationError
        Если `t` <= -273.15.
    """
    if t <= _ABSOLUTE_ZERO_C:
        raise InputValidationError(
            f'{name} должна быть выше абсолютного нуля ({_ABSOLUTE_ZERO_C} °C). Передано: {t}')


def validate_composition_normalized(composition: dict, tol: float = 1e-6) -> None:
    """
    Проверяет, что сумма мольных долей состава равна 1 (в пределах `tol`).

    Parameters
    ----------
    composition : dict
        `{имя_компонента: мольная_доля}`.
    tol : float
        Допуск на отклонение суммы от 1.

    Raises
    ------
    InputValidationError
        Если сумма долей отклоняется от 1 больше чем на `tol`.
    """
    total = sum(composition.values())
    if abs(total - 1.0) > tol:
        raise InputValidationError(
            f'Сумма мольных долей состава должна быть равна 1 (сейчас {total}). '
            f'Вызовите composition.normalize_composition() перед созданием модели.')
