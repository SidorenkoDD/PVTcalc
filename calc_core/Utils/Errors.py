"""Иерархия исключений проекта (все — плоские подклассы Exception, без своих атрибутов)."""


class NoComponentError(Exception):
    """Компонент отсутствует в БД (`Utils/DB.json`) — см. `Composition.check_and_sort_composition`."""
    pass

class LenthMissMatchError(Exception):
    """Не совпадают длины двух массивов (опечатка в имени — Length; используется в живом коде `SeparatorTest`)."""
    pass

class nStagesError(Exception):
    """Заготовлено под ошибку количества ступеней эксперимента; в текущем живом коде не используется."""
    pass

class CompositionSumError(Exception):
    """Заготовлено под ошибку суммы мольных долей состава; в текущем живом коде не используется."""
    pass

class InvalidMolarFractionError(Exception):
    """Мольная доля компонента <= 0 — см. `Composition._normalize_composition`."""
    pass

class InvalidPressureSequence(Exception):
    """Заготовлено под ошибку немонотонной последовательности давлений; в текущем живом коде не используется."""
    pass

class InvalidExcelComponentType(Exception):
    """Некорректный тип данных для компонента при загрузке состава из Excel (`Utils/CompositionLoader.py`)."""
    pass

class InvalidExcelValueType(Exception):
    """Некорректный тип значения при загрузке состава из Excel (`Utils/CompositionLoader.py`)."""
    pass

class StopIterationError(Exception):
    """Итеративный расчёт (тест стабильности, флэш) не сошёлся за отведённый лимит итераций."""
    pass

class ConvergenceError(Exception):
    """Заготовлено под общую ошибку несходимости; в текущем живом коде не используется (используется StopIterationError)."""
    pass