"""
Сохранение/восстановление состояния сессии GUI (раскладка + выбор).

Отдельный файл `session.json`, **не** `models.json` (тот — про сами модели
и результаты движка). В Фазе 0+1 хранит минимум: последнюю открытую модель
и геометрию окна. Раскладка докинга DearPyGui сохраняется штатным
механизмом самой библиотеки в отдельный ini-файл (см. `gui.view.app`); сюда
кладётся высокоуровневый выбор, не зависящий от фреймворка представления.
"""

import json
import logging
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

DEFAULT_SESSION_PATH = "gui_session.json"


@dataclass
class SessionState:
    """Фреймворк-независимый снимок сессии."""

    active_model_id: Optional[str] = None
    window_width: int = 1280
    window_height: int = 800


def load_session(path: str = DEFAULT_SESSION_PATH) -> SessionState:
    """Читает сессию; при отсутствии/повреждении файла — значения по умолчанию."""
    p = Path(path)
    if not p.exists():
        return SessionState()
    try:
        with open(p, "r", encoding="utf-8") as f:
            data = json.load(f)
        return SessionState(**{k: data[k] for k in data if k in SessionState.__annotations__})
    except (json.JSONDecodeError, TypeError) as exc:
        logger.warning("Сессия %s повреждена (%s), беру значения по умолчанию", p, exc)
        return SessionState()


def save_session(state: SessionState, path: str = DEFAULT_SESSION_PATH) -> None:
    """Пишет сессию в JSON."""
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(asdict(state), f, indent=2, ensure_ascii=False)
    logger.debug("Сессия сохранена в %s", p)
