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
from datetime import datetime
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

DEFAULT_SESSION_PATH = "gui_session.json"


@dataclass
class SessionState:
    """
    Фреймворк-независимый снимок сессии («продолжить с того же места»).

    v2: рабочие пространства хранятся **per-model** в `workspaces`
    (`{model_id: {"open_tabs", "active_tab", "flashes", "experiments"}}`) —
    каждая модель помнит свои вкладки и расчёты; восстановление происходит
    при входе в модель со страницы Projects. `active_model_id` — последняя
    активная модель (для строки «Continue last»).

    Поля `open_tabs`/`active_tab`/`flashes`/`experiments` — legacy v1
    (workspace одной активной модели); при загрузке мигрируются в
    `workspaces[active_model_id]` (см. `load_session`).
    """

    version: int = 2
    active_model_id: Optional[str] = None
    window_width: int = 1280
    window_height: int = 800
    # v2: {model_id: workspace-dict}
    workspaces: Optional[dict] = None
    # --- legacy v1 (не использовать напрямую; мигрируются при load) ---
    open_tabs: Optional[list] = None
    active_tab: Optional[str] = None
    flashes: Optional[list] = None
    experiments: Optional[list] = None
    # ISO-время последнего сохранения (проставляется в save_session)
    saved_at: Optional[str] = None


def load_session(path: str = DEFAULT_SESSION_PATH) -> SessionState:
    """
    Читает сессию; при отсутствии/повреждении файла — значения по умолчанию.
    Прозрачно мигрирует v1 → v2 (workspace одной модели → `workspaces`).
    """
    p = Path(path)
    if not p.exists():
        return SessionState()
    try:
        with open(p, "r", encoding="utf-8") as f:
            data = json.load(f)
        state = SessionState(**{k: data[k] for k in data
                                if k in SessionState.__annotations__})
    except (json.JSONDecodeError, TypeError) as exc:
        logger.warning("Сессия %s повреждена (%s), беру значения по умолчанию", p, exc)
        return SessionState()

    # миграция v1 → v2: одиночный workspace заворачивается в workspaces
    if state.workspaces is None and state.active_model_id and (
            state.flashes or state.experiments or state.open_tabs):
        state.workspaces = {state.active_model_id: {
            "open_tabs": state.open_tabs or [],
            "active_tab": state.active_tab,
            "flashes": state.flashes or [],
            "experiments": state.experiments or [],
        }}
        logger.info("Сессия мигрирована v1 -> v2 (модель %s)", state.active_model_id)
    state.open_tabs = None
    state.active_tab = None
    state.flashes = None
    state.experiments = None
    state.version = 2
    return state


def save_session(state: SessionState, path: str = DEFAULT_SESSION_PATH) -> None:
    """Пишет сессию в JSON, проставляя время сохранения (`saved_at`)."""
    state.saved_at = datetime.now().isoformat(timespec="seconds")
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "w", encoding="utf-8") as f:
        json.dump(asdict(state), f, indent=2, ensure_ascii=False)
    logger.debug("Сессия сохранена в %s", p)
