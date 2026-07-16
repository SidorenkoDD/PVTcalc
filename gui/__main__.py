"""
Точка входа GUI: `python -m gui`.

Собирает фреймворк-независимый слой (репозиторий + состояние + сессия) и
передаёт его во View на DearPyGui. Логирование настраивается здесь один раз
(библиотечный код `calc_core`/`gui` только пишет в свои логгеры).
"""

import argparse
import logging

from calc_core.Utils.Logging import configure_logging

from gui.app_state import AppState
from gui.services.model_repository import ModelRepository
from gui.session import load_session
from gui.view.app import PVTcalcApp


def main() -> None:
    parser = argparse.ArgumentParser(description="PVTcalc GUI (DearPyGui)")
    parser.add_argument("--db", default="models.json",
                        help="Путь к файлу моделей (по умолчанию models.json)")
    parser.add_argument("--log-level", default="INFO",
                        help="Уровень логирования (DEBUG/INFO/WARNING/...)")
    args = parser.parse_args()

    configure_logging(level=getattr(logging, args.log_level.upper(), logging.INFO))

    repository = ModelRepository(db_path=args.db)
    state = AppState(repository)
    session = load_session()

    app = PVTcalcApp(state, session)
    app.run()


if __name__ == "__main__":
    main()
