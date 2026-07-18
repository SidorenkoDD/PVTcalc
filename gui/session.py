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

from calc_core.Utils.AtomicFile import atomic_write_json

logger = logging.getLogger(__name__)

DEFAULT_SESSION_PATH = "gui_session.json"


@dataclass
class SessionState:
    """
    Фреймворк-независимый снимок сессии («продолжить с того же места»).

    v3: рабочие пространства хранятся **per-model** в `workspaces`; каждый
    workspace содержит общий список сериализованных узлов, открытые вкладки,
    активную вкладку и состояние сравнения. Формат v2 также читается —
    каждая модель помнит свои вкладки и расчёты; восстановление происходит
    при входе в модель со страницы Projects. `active_project_id` — последний
    проект, а `active_model_id` — последняя активная модель внутри него.

    Поля `open_tabs`/`active_tab`/`flashes`/`experiments` — legacy v1
    (workspace одной активной модели); при загрузке мигрируются в
    `workspaces[active_model_id]` (см. `load_session`).
    """

    version: int = 3
    active_project_id: Optional[str] = None
    active_model_id: Optional[str] = None
    window_width: int = 1280
    window_height: int = 800
    # v3: {model_id: workspace-dict}
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
    Прозрачно мигрирует v1 в per-model структуру; codec представления затем
    читает как старые v2-коллекции, так и общий список узлов v3.
    """
    p = Path(path)
    if not p.exists():
        return SessionState()
    try:
        with open(p, "r", encoding="utf-8") as f:
            data = json.load(f)
        if not isinstance(data, dict):
            raise TypeError("корень сессии должен быть JSON-объектом")

        def optional_str(value):
            return value if isinstance(value, str) else None

        def window_size(value, default: int) -> int:
            if isinstance(value, bool) or not isinstance(value, (int, float)):
                return default
            return max(320, min(10000, int(value)))

        raw_workspaces = data.get("workspaces")
        workspaces = None
        if isinstance(raw_workspaces, dict):
            workspaces = {
                key: value for key, value in raw_workspaces.items()
                if isinstance(key, str) and isinstance(value, dict)
            }
        state = SessionState(
            version=3,
            active_project_id=optional_str(data.get("active_project_id")),
            active_model_id=optional_str(data.get("active_model_id")),
            window_width=window_size(data.get("window_width"), 1280),
            window_height=window_size(data.get("window_height"), 800),
            workspaces=workspaces,
            open_tabs=(data.get("open_tabs")
                       if isinstance(data.get("open_tabs"), list) else None),
            active_tab=optional_str(data.get("active_tab")),
            flashes=(data.get("flashes")
                     if isinstance(data.get("flashes"), list) else None),
            experiments=(data.get("experiments")
                         if isinstance(data.get("experiments"), list) else None),
            saved_at=optional_str(data.get("saved_at")),
        )
    except (json.JSONDecodeError, TypeError, ValueError, OverflowError) as exc:
        logger.warning("Сессия %s повреждена (%s), беру значения по умолчанию", p, exc)
        return SessionState()

    # Миграция v1: одиночный workspace заворачивается в per-model workspaces.
    if state.workspaces is None and state.active_model_id and (
            state.flashes or state.experiments or state.open_tabs):
        state.workspaces = {state.active_model_id: {
            "open_tabs": state.open_tabs or [],
            "active_tab": state.active_tab,
            "flashes": state.flashes or [],
            "experiments": state.experiments or [],
        }}
        logger.info("Сессия мигрирована v1 -> v3 (модель %s)", state.active_model_id)
    state.open_tabs = None
    state.active_tab = None
    state.flashes = None
    state.experiments = None
    state.version = 3
    return state


def save_session(state: SessionState, path: str = DEFAULT_SESSION_PATH) -> None:
    """Атомарно пишет сессию в JSON, проставляя время сохранения."""
    state.saved_at = datetime.now().isoformat(timespec="seconds")
    p = Path(path)
    atomic_write_json(p, asdict(state))
    logger.debug("Сессия сохранена в %s", p)
