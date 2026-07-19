"""
Главное окно GUI на DearPyGui (тонкий View над `gui.app_state.AppState`).

Навигация в стиле IDE:
- слева — раскрывающееся дерево `Project → Model → (Composition / Flash-runs /
  …)`; текущий узел подсвечен; состояние разворота узлов хранится во View
  (`_expanded_models`/`_expanded_cats`);
- центральная рабочая область — **набор вкладок** (как редакторы файлов):
  клик по узлу в дереве открывает/фокусирует его вкладку, крестик закрывает.
  Composition открывается вкладкой (с внутренними под-вкладками Properties/BIP),
  каждый флэш — своей вкладкой; история всех флэшей видна в дереве.

Весь видимый пользователю текст — на английском; комментарии/докстринги — на
русском (конвенция проекта). Тема — светлая (`gui/view/theme.py`).
"""

import logging

import dearpygui.dearpygui as dpg

from gui.app_state import AppState, NodeKind, StateChange, StateChangeKind
from gui.calculation_coordinator import CalculationCoordinator
from gui.session import SessionState, save_session
from gui.view.composition_view import CompositionViewMixin
from gui.view.contracts import ViewContext
from gui.view.dialogs_view import DialogsViewMixin
from gui.view.envelope_view import EnvelopeViewMixin
from gui.view.experiment_view import ExperimentViewMixin
from gui.view.flash_view import FlashViewMixin
from gui.view.new_fluid_form import NewFluidForm
from gui.view.projects_view import ProjectsViewMixin
from gui.view.theme import build_light_theme
from gui.view.workspace_view import WorkspaceViewMixin
from gui.workspace_codec import restore_workspace, snapshot_workspace

logger = logging.getLogger(__name__)

# Теги контейнеров
_PRIMARY = "primary_window"
_TREE_PANEL = "tree_panel"
_WORKSPACE_PANEL = "workspace_panel"
_MODEL_TREE = "model_tree"
_WORKSPACE = "workspace_content"
_STATUS_BAR = "status_bar"
# статичные контейнеры-экраны (создаются один раз, переключаются show=)
_PROJECTS_SCREEN = "projects_screen"
_PROJECTS_CONTENT = "projects_content"
_WORKSPACE_SCREEN = "workspace_screen"

_STATUS_H = 28    # высота строки статуса
_TREE_W = 320     # ширина левой панели

class PVTcalcApp(
    ProjectsViewMixin,
    CompositionViewMixin,
    FlashViewMixin,
    ExperimentViewMixin,
    EnvelopeViewMixin,
    DialogsViewMixin,
    WorkspaceViewMixin,
):
    """View-контроллер: связывает `AppState` с виджетами DearPyGui."""

    def __init__(self, state: AppState, session: SessionState):
        self._view_context = ViewContext(
            state=state,
            session=session,
            jobs=CalculationCoordinator(),
        )
        # состояние разворота дерева (во View, чтобы переживать перерисовку)
        self._expanded_models: set[str] = set()
        self._expanded_cats: set[str] = set()
        # выбор узлов для сравнения (Ctrl+клик в дереве)
        # Полный адрес флэша: локальные node_id (например, flash_1) могут
        # повторяться в разных моделях одного проекта.
        self._compare_selection: set[tuple[str, str]] = set()
        # id таб-бара и вкладок рабочей области (захватываем при отрисовке,
        # чтобы не плодить фиксированные алиасы при пересоздании)
        self._tabbar_id: int | None = None
        self._tab_ids: dict[str, int] = {}
        self._tab_content_ids: dict[str, int] = {}
        self._workspace_crumb_id: int | None = None
        # id полей ввода (захватываем при отрисовке, без фиксированных алиасов)
        self._bip_ids: dict[tuple[int, int], int] = {}
        self._flash_input_ids: dict[str, tuple[int, int]] = {}
        self._exp_input_ids: dict[str, dict] = {}
        self._exp_chart_holder: dict[str, int] = {}  # контейнер графиков вкладки
        self._lab_data_holder: dict[str, int] = {}
        self._lab_data_controls: dict[str, tuple[int, int, int]] = {}
        # модальное окно параметров фазовой огибающей (собирается перед расчётом)
        self._env_dialog_win: int | None = None
        self._env_dialog_ids: dict = {}
        self._env_dialog_node: str | None = None  # None = новый узел, иначе re-run
        # окно экспорта модели во внешний формат (E300 и т.п.)
        self._export_win: int | None = None
        self._export_ids: dict = {}
        self._export_comp = None
        self._export_label: str = "model"
        self._export_fmt: str = "e300"
        self._export_eos: str = "MPR"
        # окно настроек (константы/условия/критерии сходимости)
        self._settings_win: int | None = None
        self._settings_ids: dict = {}
        # debounce автосохранения модели/сессии (номер последнего поколения)
        self._model_save_generation: int = 0
        self._session_save_generation: int = 0
        # модели, чей workspace уже восстановлен из сессии в этом запуске
        self._restored_models: set[str] = set()
        # стек открытых модальных окон (для закрытия по Esc, верхнее первым)
        self._modals: list[int] = []
        # экран Projects: выбранный проект + распознавание двойного клика
        self._selected_project: str | None = None
        self._proj_row_ids: dict[str, int] = {}
        self._last_proj_click: str | None = None
        self._last_proj_click_time: float = 0.0
        # импорт Excel: путь, id окна предпросмотра и его виджетов
        self._excel_path: str | None = None
        self._excel_win: int | None = None
        self._excel_sheet_id: int | None = None
        self._excel_header_id: int | None = None
        self._excel_preview_group: int | None = None
        # импорт E300: обязательный предпросмотр до изменения базы
        self._e300_path: str | None = None
        self._e300_win: int | None = None
        # форма создания нового флюида — модальное окно
        self._new_fluid_win: int | None = None
        # модальное окно копирования модели из рабочего дерева
        self._duplicate_win: int | None = None
        self._duplicate_ids: dict[str, int] = {}
        self._new_fluid_form = NewFluidForm(
            get_db_path=lambda: self._state.db_path,
            on_created=self._on_fluid_created,
            on_cancel=self._close_new_fluid_modal,
            set_status=self._set_status,
        )
        state.subscribe_changes(self._render_change)
        state.subscribe(self._schedule_session_autosave)
        state.subscribe_dirty(self._schedule_model_autosave)

    # --- жизненный цикл --------------------------------------------------

    def run(self) -> None:
        dpg.create_context()
        dpg.bind_theme(build_light_theme())

        dpg.create_viewport(title="PVTcalc", width=self._session.window_width,
                            height=self._session.window_height)
        self._build_menu()
        self._build_layout()
        self._build_shortcuts()

        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.set_primary_window(_PRIMARY, True)

        self._state.refresh_model_list()  # notify -> рендер экрана Projects

        dpg.start_dearpygui()

        self._save_dirty_models()
        self._persist_session()
        dpg.destroy_context()

    # --- статичный каркас ------------------------------------------------

    def _build_menu(self) -> None:
        with dpg.viewport_menu_bar():
            with dpg.menu(label="Project"):
                dpg.add_menu_item(label="Save model   (Ctrl+S)",
                                  callback=lambda: self._on_save_model())
                dpg.add_menu_item(label="Projects home",
                                  callback=self._on_back_to_projects)
                dpg.add_menu_item(label="Refresh models",
                                  callback=lambda: self._state.refresh_model_list())
            with dpg.menu(label="Edit"):
                dpg.add_menu_item(label="Undo   (Ctrl+Z)",
                                  callback=lambda: self._do_undo())
                dpg.add_menu_item(label="Redo   (Ctrl+Y)",
                                  callback=lambda: self._do_redo())
            with dpg.menu(label="Export"):
                dpg.add_menu_item(label="Export model...",
                                  callback=self._on_open_export_dialog)
            with dpg.menu(label="Settings"):
                dpg.add_menu_item(label="Engine constants (read only)...",
                                  callback=self._on_open_settings)

    def _build_shortcuts(self) -> None:
        """Глобальные горячие клавиши: Del (в дереве), Ctrl+Z / Ctrl+Y."""
        with dpg.handler_registry():
            dpg.add_key_press_handler(dpg.mvKey_Delete, callback=self._on_key_delete)
            dpg.add_key_press_handler(dpg.mvKey_Z, callback=self._on_key_z)
            dpg.add_key_press_handler(dpg.mvKey_Y, callback=self._on_key_y)
            dpg.add_key_press_handler(dpg.mvKey_S, callback=self._on_key_s)
            dpg.add_key_press_handler(dpg.mvKey_Return, callback=self._on_key_enter)
            npad = getattr(dpg, "mvKey_NumPadEnter", None) or getattr(
                dpg, "mvKey_KeyPadEnter", None)
            if npad is not None:
                dpg.add_key_press_handler(npad, callback=self._on_key_enter)
            dpg.add_key_press_handler(dpg.mvKey_Escape, callback=self._on_key_escape)
            dpg.add_key_press_handler(dpg.mvKey_Up, callback=self._on_key_up)
            dpg.add_key_press_handler(dpg.mvKey_Down, callback=self._on_key_down)

    # --- модальные окна (закрытие по Esc) -------------------------------

    def _track_modal(self, win: int) -> int:
        """Регистрирует модальное окно в стеке для закрытия по Esc. Возвращает id."""
        self._modals.append(win)
        return win

    def _has_open_modal(self) -> bool:
        return any(dpg.does_item_exist(w) for w in self._modals)

    def _on_key_escape(self, sender, app_data) -> None:
        """Esc закрывает верхнее открытое модальное окно."""
        while self._modals:
            win = self._modals[-1]
            if dpg.does_item_exist(win):
                dpg.delete_item(win)
                self._modals.pop()
                return
            self._modals.pop()

    def _on_key_up(self, sender, app_data) -> None:
        self._project_move_selection(-1)

    def _on_key_down(self, sender, app_data) -> None:
        self._project_move_selection(+1)

    def _project_move_selection(self, delta: int) -> None:
        """Стрелки ↑/↓ двигают выбор в таблице проектов (когда нет модалок)."""
        if self._state.current_screen != "projects" or self._has_open_modal():
            return
        ids = list(self._proj_row_ids.keys())
        if not ids:
            return
        if self._selected_project in ids:
            i = min(len(ids) - 1, max(0, ids.index(self._selected_project) + delta))
        else:
            i = 0
        mid = ids[i]
        self._selected_project = mid
        for m, sid in self._proj_row_ids.items():
            if dpg.does_item_exist(sid):
                dpg.set_value(sid, m == mid)
        self._set_status(f"Selected '{mid}' (double-click or Enter to open).")

    @staticmethod
    def _ctrl_down() -> bool:
        return (dpg.is_key_down(dpg.mvKey_LControl)
                or dpg.is_key_down(dpg.mvKey_RControl))

    @staticmethod
    def _shift_down() -> bool:
        return (dpg.is_key_down(dpg.mvKey_LShift)
                or dpg.is_key_down(dpg.mvKey_RShift))

    def _on_key_delete(self, sender, app_data) -> None:
        if self._state.current_screen != "workspace":
            return
        # только когда курсор над панелью дерева (не мешаем вводу в поля)
        if not (dpg.is_item_hovered(_TREE_PANEL) or dpg.is_item_hovered(_MODEL_TREE)):
            return
        node = self._state.active_node
        # удаляем любой узел-расчёт (Composition не удаляется — см. delete_node)
        if node is not None and node.kind is not NodeKind.COMPOSITION:
            self._state.delete_node(node.node_id)
            self._set_status(f"{node.kind.name.capitalize()} node deleted (Del).")

    def _on_key_z(self, sender, app_data) -> None:
        if self._state.current_screen != "workspace":
            return
        if self._ctrl_down():
            self._do_redo() if self._shift_down() else self._do_undo()

    def _on_key_y(self, sender, app_data) -> None:
        if self._state.current_screen != "workspace":
            return
        if self._ctrl_down():
            self._do_redo()

    def _on_key_s(self, sender, app_data) -> None:
        if self._ctrl_down():
            self._on_save_model()

    def _on_key_enter(self, sender, app_data) -> None:
        # Enter на странице Projects открывает выбранную модель (не при модалке)
        if self._state.current_screen != "projects" or self._has_open_modal():
            return
        project_id = self._selected_project
        if project_id and project_id in self._state.projects:
            self._open_project(project_id)

    def _do_undo(self) -> None:
        if self._state.can_undo():
            self._state.undo()
            self._set_status("Undo.")
        else:
            self._set_status("Nothing to undo.")

    def _do_redo(self) -> None:
        if self._state.can_redo():
            self._state.redo()
            self._set_status("Redo.")
        else:
            self._set_status("Nothing to redo.")

    def _on_save_model(self) -> None:
        model = self._state.active_model
        if model is None or not model.loaded:
            self._set_status("No loaded model to save.")
            return
        try:
            self._state.save_model(model.model_id)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Не удалось сохранить модель %s", model.model_id)
            self._set_status(f"Save failed: {exc}")
            return
        self._set_status(f"Model '{model.model_id}' saved.")

    def _save_dirty_models(self) -> None:
        """Сохраняет все dirty-модели; используется debounce и shutdown."""
        try:
            saved = self._state.save_all_dirty()
        except Exception as exc:  # noqa: BLE001
            logger.exception("Автосохранение моделей не удалось")
            self._set_status(f"Autosave failed: {exc}")
            return
        if saved:
            self._set_status("Autosaved: " + ", ".join(saved))

    def _schedule_model_autosave(self, model_id: str) -> None:
        """Debounce сохранения правок состава примерно на одну секунду."""
        self._model_save_generation += 1
        generation = self._model_save_generation

        def save_if_latest(*_args) -> None:
            if generation == self._model_save_generation:
                self._save_dirty_models()

        dpg.set_frame_callback(dpg.get_frame_count() + 60, save_if_latest)

    def _schedule_session_autosave(self) -> None:
        """Debounce снимка UI-сессии после изменений состояния."""
        self._session_save_generation += 1
        generation = self._session_save_generation

        def save_if_latest(*_args) -> None:
            if generation == self._session_save_generation:
                self._persist_session()

        dpg.set_frame_callback(dpg.get_frame_count() + 60, save_if_latest)

    def _build_layout(self) -> None:
        """
        Три экрана-контейнера (Projects / New fluid / Workspace), созданные
        один раз и переключаемые `show=` (см. `_render`); глобальная строка
        статуса внизу. Внутри экранов — прежний паттерн перерисовки
        (delete children + пересоздание с захватом id).
        """
        with dpg.window(tag=_PRIMARY, no_title_bar=True, no_move=True,
                        no_resize=True, no_collapse=True):
            # --- экран Projects (стартовый) ---
            with dpg.child_window(tag=_PROJECTS_SCREEN, height=-_STATUS_H,
                                  border=False, show=True):
                dpg.add_group(tag=_PROJECTS_CONTENT)

            # --- экран рабочего пространства (дерево + вкладки) ---
            with dpg.group(tag=_WORKSPACE_SCREEN, show=False):
                # таблица-сплиттер: перетаскиваемая граница между деревом и
                # рабочей областью (ширина левой панели меняется мышью)
                with dpg.table(header_row=False, resizable=True,
                               borders_innerV=True,
                               policy=dpg.mvTable_SizingStretchProp,
                               height=-_STATUS_H):
                    dpg.add_table_column(init_width_or_weight=0.24)  # дерево
                    dpg.add_table_column(init_width_or_weight=0.76)  # раб. область
                    with dpg.table_row():
                        with dpg.child_window(tag=_TREE_PANEL, width=-1,
                                              height=-1, border=True):
                            dpg.add_button(label="< Projects",
                                           callback=self._on_back_to_projects)
                            dpg.add_separator()
                            dpg.add_group(tag=_MODEL_TREE)
                        with dpg.child_window(tag=_WORKSPACE_PANEL, width=-1,
                                              height=-1, border=True):
                            dpg.add_group(tag=_WORKSPACE)

            dpg.add_separator()
            dpg.add_text("Ready.", tag=_STATUS_BAR)

    # --- полная перерисовка ---------------------------------------------

    def _render(self) -> None:
        """Диспетчер по текущему экрану: показывает нужный контейнер и
        перерисовывает только его содержимое."""
        screen = self._state.current_screen
        if dpg.does_item_exist(_PROJECTS_SCREEN):
            dpg.configure_item(_PROJECTS_SCREEN, show=(screen == "projects"))
            dpg.configure_item(_WORKSPACE_SCREEN, show=(screen == "workspace"))
        if screen == "projects":
            self._render_projects()
        elif screen == "workspace":
            self._render_tree()
            self._render_workspace()
        # Форма нового флюида — модальное окно поверх Projects (рендерится один
        # раз при открытии, не по notify: держит ввод локально, не теряет фокус).

    def _render_change(self, change: StateChange) -> None:
        """Обновляет только затронутую область после команды AppState.

        Полная пересборка остаётся явным fallback для публичного ``notify()``.
        Особенно важен путь ``NODE``: завершение фонового расчёта не трогает
        соседние вкладки и не сбрасывает их локальные DPG-состояния.
        """
        if change.kind is StateChangeKind.FULL:
            self._render()
            return
        if not dpg.does_item_exist(_PROJECTS_SCREEN):
            return

        screen = self._state.current_screen
        dpg.configure_item(_PROJECTS_SCREEN, show=(screen == "projects"))
        dpg.configure_item(_WORKSPACE_SCREEN, show=(screen == "workspace"))

        if screen == "projects":
            if change.kind in (StateChangeKind.NAVIGATION,
                               StateChangeKind.MODEL_LIST):
                self._render_projects()
            return

        if change.kind is StateChangeKind.NAVIGATION:
            self._render_tree()
            self._render_workspace()
        elif change.kind is StateChangeKind.MODEL_LIST:
            self._render_tree()
            if self._state.active_model is None:
                self._render_workspace()
        elif change.kind is StateChangeKind.TREE:
            self._render_tree()
        elif change.kind is StateChangeKind.WORKSPACE:
            self._render_tree()
            self._render_workspace()
        elif change.kind is StateChangeKind.NODE:
            ref = change.node_ref
            if ref is None or ref.model_id == self._state.active_model_id:
                self._render_tree()
                variant = self._state.active_variant
                is_open = (ref is not None and variant is not None
                           and ref.node_id in variant.open_node_ids)
                if ref is None:
                    self._render_workspace()
                elif is_open and not self._render_node_content(ref.node_id):
                    self._render_workspace()

    # ==================================================================
    #  Экран Projects (стартовый)
    # ==================================================================
    @staticmethod
    def _fmt(value) -> str:
        if value is None:
            return "-"
        if isinstance(value, (int, float)):
            return f"{value:.5g}"
        return str(value)

    @staticmethod
    def _g(value) -> str:
        return f"{value:g}" if isinstance(value, (int, float)) else str(value)

    def _theme_stale(self):
        """Ленивая тема (тёмно-оранжевый текст) для пометки устаревшего/ошибки."""
        if getattr(self, "_stale_theme_id", None) is None:
            with dpg.theme() as theme_id:
                with dpg.theme_component(dpg.mvText):
                    dpg.add_theme_color(dpg.mvThemeCol_Text, (200, 110, 0))
            self._stale_theme_id = theme_id
        return self._stale_theme_id

    def _set_status(self, text: str) -> None:
        if dpg.does_item_exist(_STATUS_BAR):
            dpg.set_value(_STATUS_BAR, text)

    # ==================================================================
    #  Сессия
    # ==================================================================

    @staticmethod
    def _format_saved_at(iso: str | None) -> str:
        """ISO-строку `saved_at` -> 'YYYY-MM-DD HH:MM' (для строки Continue last)."""
        if not iso:
            return ""
        try:
            from datetime import datetime
            return datetime.fromisoformat(iso).strftime("%Y-%m-%d %H:%M")
        except (ValueError, TypeError):
            return ""

    def _restore_workspace(self, model_id: str) -> None:
        """
        Лениво восстанавливает workspace модели из `session.workspaces[mid]`
        (узлы/результаты/вкладки/сравнение). Повторный вызов для той же модели —
        no-op (guard `_restored_models`). Модель должна быть активной.
        """
        if model_id in self._restored_models:
            return
        self._restored_models.add(model_id)
        ws = (self._session.workspaces or {}).get(model_id) or {}
        variant = self._state.active_variant
        if variant is None or not ws:
            return

        restore_workspace(variant, ws)

        # Развернуть только реально заполненные категории восстановленного дерева.
        self._expanded_models.add(model_id)
        if variant.flash_runs():
            self._expanded_cats.add(f"{model_id}:flash")
        if variant.experiment_runs():
            self._expanded_cats.add(f"{model_id}:exp")
            for node in variant.experiment_runs():
                kind = str(node.params.get("kind", "")).lower()
                if kind:
                    self._expanded_cats.add(f"{model_id}:exp:{kind}")
        if variant.envelope_runs():
            self._expanded_cats.add(f"{model_id}:env")

        self._state.clear_history()  # восстановление не откатываем

    @staticmethod
    def _workspace_snapshot(variant) -> dict:
        """Совместимый wrapper над framework-independent codec workspace."""
        return snapshot_workspace(variant)

    def _persist_session(self) -> None:
        self._session.active_project_id = self._state.active_project_id
        self._session.active_model_id = self._state.active_model_id
        try:
            self._session.window_width = dpg.get_viewport_width()
            self._session.window_height = dpg.get_viewport_height()
        except Exception:  # noqa: BLE001
            pass

        # merge: обновляем workspaces загруженных моделей, не трогая остальные
        if self._session.workspaces is None:
            self._session.workspaces = {}
        for mid, model in self._state.models.items():
            variant = model.variants.get("base") if model.loaded else None
            if variant is not None:
                self._session.workspaces[mid] = self._workspace_snapshot(variant)

        save_session(self._session)
