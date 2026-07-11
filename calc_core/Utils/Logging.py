"""Единая точка настройки логирования для PVTcalc.

Библиотечный код (calc_core/**) сам никогда не вызывает logging.basicConfig()
и не навешивает handlers — модуль только объявляет
`logger = logging.getLogger(__name__)` и пишет через него. `__name__`
внутри пакета даёт иерархию вида `calc_core.VLE.Flash`, `calc_core.EOS.BrusilovskiyEOS`
и т.д., корень этой иерархии — логгер с именем `calc_core`.

Настройка (куда писать, какой уровень) — прерогатива того, кто использует
движок: вызвать configure_logging() один раз в начале ноутбука/скрипта/теста,
а не полагаться на побочный эффект чьего-то импорта.
"""

import logging

_ENGINE_LOGGER_NAME = "calc_core"


def configure_logging(level: int = logging.INFO, fmt: str | None = None) -> logging.Logger:
    """Настраивает вывод логов всего движка (calc_core.*) в консоль.

    Идемпотентна — повторный вызов не создаёт дублирующихся handlers.
    """
    logger = logging.getLogger(_ENGINE_LOGGER_NAME)
    logger.setLevel(level)

    if not logger.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(
            logging.Formatter(
                fmt or "%(asctime)s | %(name)s | %(levelname)s | %(message)s",
                datefmt="%Y-%m-%d %H:%M:%S",
            )
        )
        logger.addHandler(handler)

    return logger
