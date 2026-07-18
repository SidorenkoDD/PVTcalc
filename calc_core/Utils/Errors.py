"""
Иерархия исключений проекта.

Все исключения, которые движок бросает осознанно, наследуются от базового
`PVTCalcError` — чтобы внешний код (в первую очередь UI/сервисный слой) мог
одним `except PVTCalcError` отделить «доменную» ошибку движка от неожиданного
бага (обычного `Exception`). Две смысловые ветки:

- `InputValidationError` — некорректный вход от вызывающего кода/пользователя
  (то, что UI показывает пользователю: «поправьте состав/давление/температуру»).
- `ConvergenceError` — численный метод не сошёлся (то, что UI показывает как
  «расчёт не сошёлся», а не как некорректный ввод).

Имена классов и пути импорта сохранены прежними — переподчинение меняет только
базовый класс, поэтому существующий `except NoComponentError`/`except
StopIterationError` продолжает работать, а `except PVTCalcError` теперь ловит
их все.
"""


class PVTCalcError(Exception):
    """Базовое исключение движка — всё, что PVTcalc бросает осознанно."""
    pass


class InputValidationError(PVTCalcError):
    """Некорректный вход от вызывающего кода/пользователя (см. `Utils/Validation.py`)."""
    pass


class ConvergenceError(PVTCalcError):
    """Численный метод не сошёлся (поднимается фасадами по сигналу калькулятора, см. `PhaseEnvelopeFacade`)."""
    pass


class NoComponentError(InputValidationError):
    """Компонент отсутствует в БД (`Utils/DB.json`) — см. `Composition.check_and_sort_composition`."""
    pass


class InvalidMolarFractionError(InputValidationError):
    """Мольная доля компонента <= 0 — см. `Composition._normalize_composition`."""
    pass


class CompositionSumError(InputValidationError):
    """Сумма мольных долей состава != 1 — см. `Utils/Validation.validate_composition_normalized`."""
    pass


class StopIterationError(ConvergenceError):
    """Итеративный расчёт (тест стабильности, флэш) не сошёлся за отведённый лимит итераций."""
    pass


class LengthMismatchError(PVTCalcError):
    """Не совпадают длины двух массивов."""
    pass


# Обратная совместимость с исторической опечаткой.
LenthMissMatchError = LengthMismatchError


class nStagesError(PVTCalcError):
    """Заготовлено под ошибку количества ступеней эксперимента; в текущем живом коде не используется."""
    pass


class InvalidPressureSequence(PVTCalcError):
    """Заготовлено под ошибку немонотонной последовательности давлений; в текущем живом коде не используется."""
    pass


class InvalidExcelComponentType(PVTCalcError):
    """Некорректный тип данных для компонента при загрузке состава из Excel (`Utils/CompositionLoader.py`)."""
    pass


class InvalidExcelValueType(PVTCalcError):
    """Некорректный тип значения при загрузке состава из Excel (`Utils/CompositionLoader.py`)."""
    pass
