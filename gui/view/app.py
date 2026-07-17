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
import threading
import time
from pathlib import Path

import dearpygui.dearpygui as dpg

from gui.app_state import AppState, NodeKind, NodeStatus
from gui.services import composition_service as comp_svc
from gui.services import experiment_service as exp_svc
from gui.services import flash_service
from gui.services import project_service as proj_svc
from gui.session import SessionState, save_session
from gui.view.new_fluid_form import NewFluidForm
from gui.view.theme import build_light_theme

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
_BIP_CELL_W = 60  # ширина ячейки матрицы BIP

# Скалярные свойства компонент в таблице состава (ключ -> заголовок колонки).
_COMPOSITION_COLUMNS: list[tuple[str, str]] = [
    ("molar_mass", "M, g/mol"),
    ("critical_temperature", "Tc, K"),
    ("critical_pressure", "Pc, bar"),
    ("acentric_factor", "omega"),
    ("critical_volume", "Vc"),
    ("shift_parameter", "shift"),
    ("Kw", "Kw"),
]


class PVTcalcApp:
    """View-контроллер: связывает `AppState` с виджетами DearPyGui."""

    def __init__(self, state: AppState, session: SessionState):
        self._state = state
        self._session = session
        # состояние разворота дерева (во View, чтобы переживать перерисовку)
        self._expanded_models: set[str] = set()
        self._expanded_cats: set[str] = set()
        # выбор узлов для сравнения (Ctrl+клик в дереве)
        self._compare_selection: set[str] = set()
        # активная фоновая задача флэша (расчёт в отдельном потоке)
        self._active_job: dict | None = None
        # id таб-бара и вкладок рабочей области (захватываем при отрисовке,
        # чтобы не плодить фиксированные алиасы при пересоздании)
        self._tabbar_id: int | None = None
        self._tab_ids: dict[str, int] = {}
        # id полей ввода (захватываем при отрисовке, без фиксированных алиасов)
        self._bip_ids: dict[tuple[int, int], int] = {}
        self._flash_input_ids: dict[str, tuple[int, int]] = {}
        self._exp_input_ids: dict[str, dict] = {}
        self._exp_chart_holder: dict[str, int] = {}  # контейнер графиков вкладки
        # модели, чей workspace уже восстановлен из сессии в этом запуске
        self._restored_models: set[str] = set()
        # стек открытых модальных окон (для закрытия по Esc, верхнее первым)
        self._modals: list[int] = []
        # экран Projects: выбранная модель + распознавание двойного клика
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
        # форма создания нового флюида — модальное окно
        self._new_fluid_win: int | None = None
        self._new_fluid_form = NewFluidForm(
            get_db_path=lambda: self._state.db_path,
            on_created=self._on_fluid_created,
            on_cancel=self._close_new_fluid_modal,
            set_status=self._set_status,
        )
        state.subscribe(self._render)

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

        self._persist_session()
        dpg.destroy_context()

    # --- статичный каркас ------------------------------------------------

    def _build_menu(self) -> None:
        with dpg.viewport_menu_bar():
            with dpg.menu(label="Project"):
                dpg.add_menu_item(label="Projects home",
                                  callback=self._on_back_to_projects)
                dpg.add_menu_item(label="Refresh models",
                                  callback=lambda: self._state.refresh_model_list())
            with dpg.menu(label="Edit"):
                dpg.add_menu_item(label="Undo   (Ctrl+Z)",
                                  callback=lambda: self._do_undo())
                dpg.add_menu_item(label="Redo   (Ctrl+Y)",
                                  callback=lambda: self._do_redo())

    def _build_shortcuts(self) -> None:
        """Глобальные горячие клавиши: Del (в дереве), Ctrl+Z / Ctrl+Y."""
        with dpg.handler_registry():
            dpg.add_key_press_handler(dpg.mvKey_Delete, callback=self._on_key_delete)
            dpg.add_key_press_handler(dpg.mvKey_Z, callback=self._on_key_z)
            dpg.add_key_press_handler(dpg.mvKey_Y, callback=self._on_key_y)
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

    def _on_key_enter(self, sender, app_data) -> None:
        # Enter на странице Projects открывает выбранную модель (не при модалке)
        if self._state.current_screen != "projects" or self._has_open_modal():
            return
        mid = self._selected_project
        if mid and mid in self._state.models:
            self._open_project(mid)

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

    # ==================================================================
    #  Экран Projects (стартовый)
    # ==================================================================

    def _render_projects(self) -> None:
        if not dpg.does_item_exist(_PROJECTS_CONTENT):
            return
        dpg.delete_item(_PROJECTS_CONTENT, children_only=True)
        parent = _PROJECTS_CONTENT

        dpg.add_spacer(height=6, parent=parent)
        dpg.add_text("Projects", parent=parent)
        dpg.add_separator(parent=parent)

        # строка «продолжить последнюю сессию»
        last = self._session.active_model_id
        if last and last in self._state.models:
            saved = self._format_saved_at(self._session.saved_at)
            with dpg.group(horizontal=True, parent=parent):
                dpg.add_text(f"Continue last: {last}"
                             + (f"  (saved {saved})" if saved else ""))
                dpg.add_button(label="Open", user_data=last,
                               callback=self._on_open_model)
            dpg.add_separator(parent=parent)

        dpg.add_text("Double-click (or Enter) to open a model, right-click for "
                     "more actions.", parent=parent)

        # таблица моделей — строки выбираемые (клик = выбор, двойной = открыть)
        self._proj_row_ids = {}
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True,
                       borders_outerV=True, resizable=True, scrollY=True,
                       height=380, freeze_rows=1):
            for label in ("Model", "Field", "EOS", "Components", "T res, K",
                          "Calculated", "Created"):
                dpg.add_table_column(label=label)
            workspaces = self._session.workspaces or {}
            for model in self._state.models.values():
                s = model.summary
                with dpg.table_row():
                    sel = dpg.add_selectable(
                        label=model.title, span_columns=True,
                        default_value=(model.model_id == self._selected_project),
                        user_data=model.model_id, callback=self._on_project_click)
                    self._proj_row_ids[model.model_id] = sel
                    self._attach_project_context_menu(sel, model.model_id)
                    dpg.add_text(model.field_name or "-")
                    dpg.add_text(model.eos or "-")
                    dpg.add_text(str(model.n_components))
                    dpg.add_text(self._g(s.t_res) if s and s.t_res else "-")
                    dpg.add_text(self._calc_summary_text(model, workspaces))
                    dpg.add_text((s.created_at or "-")[:10] if s else "-")

        dpg.add_spacer(height=8, parent=parent)
        dpg.add_text("New composition", parent=parent)
        with dpg.group(horizontal=True, parent=parent):
            dpg.add_button(label="Create manually",
                           callback=self._on_new_fluid_manual)
            dpg.add_button(label="Import Excel", callback=self._on_import_excel)
            e300 = dpg.add_button(label="Import E300 (soon)", enabled=False)
            with dpg.tooltip(e300):
                dpg.add_text("Not available yet - E300 parser is planned.")

    def _calc_summary_text(self, model, workspaces: dict) -> str:
        """Текст колонки Calculated: persisted + сессионные расчёты."""
        if model.summary is None:
            return "-"
        info = proj_svc.calc_summary(model.summary,
                                     workspaces.get(model.model_id))
        parts: list[str] = []
        if info["flashes"]:
            parts.append(f"{info['flashes']} flash")
        for kind, n in info["exp_kinds"].items():
            parts.append(f"{n} {kind}")
        if info["persisted"]:
            parts.append(f"{info['persisted']} stored")
        return ", ".join(parts) if parts else "-"

    # --- навигация ---------------------------------------------------------

    def _on_open_model(self, sender, app_data, user_data) -> None:
        self._open_project(user_data)

    def _on_project_click(self, sender, app_data, user_data) -> None:
        """Клик по строке: выбор + подсветка; двойной клик — открыть модель."""
        mid = user_data
        now = time.monotonic()
        is_double = (self._last_proj_click == mid
                     and now - self._last_proj_click_time < 0.4)
        self._last_proj_click = mid
        self._last_proj_click_time = now
        self._selected_project = mid
        # подсветка без перерисовки (иначе двойной клик не долетит до строки)
        for m, sid in self._proj_row_ids.items():
            if dpg.does_item_exist(sid):
                dpg.set_value(sid, m == mid)
        if is_double:
            self._open_project(mid)
        else:
            self._set_status(f"Selected '{mid}' (double-click or Enter to open).")

    def _open_project(self, mid: str) -> None:
        if self._active_job is not None:
            self._set_status("Calculation in progress - wait or cancel first.")
            return
        try:
            self._state.enter_model(mid)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Не удалось открыть модель %s", mid)
            self._set_status(f"Failed to open '{mid}': {exc}")
            return
        self._expanded_models.add(mid)
        self._restore_workspace(mid)
        self._state._notify()
        self._set_status(f"Model '{mid}' opened.")

    # --- контекстное меню модели (просмотр состава / удаление) -------------

    def _attach_project_context_menu(self, item_id, model_id: str) -> None:
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text(model_id)
            dpg.add_separator()
            dpg.add_menu_item(label="Open", user_data=model_id,
                              callback=self._on_open_model)
            dpg.add_menu_item(label="View composition (read-only)",
                              user_data=model_id, callback=self._on_view_composition)
            dpg.add_menu_item(label="Delete model...", user_data=model_id,
                              callback=self._on_delete_model_confirm)

    def _on_view_composition(self, sender, app_data, user_data) -> None:
        mid = user_data
        try:
            comp = self._state.peek_composition(mid)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Не удалось прочитать состав %s", mid)
            self._set_status(f"Failed to read composition of '{mid}': {exc}")
            return
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        with dpg.window(label=f"Composition (read-only) - {mid}", modal=True,
                        width=720, height=560,
                        pos=(max(0, w // 2 - 360), max(0, h // 2 - 280))) as win:
            self._track_modal(win)
            data = comp.composition_data
            with dpg.table(header_row=True, borders_innerH=True, borders_outerH=True,
                           borders_innerV=True, borders_outerV=True,
                           scrollY=True, height=-40, freeze_rows=1):
                dpg.add_table_column(label="Component")
                dpg.add_table_column(label="zi, frac.")
                for _k, title in _COMPOSITION_COLUMNS:
                    dpg.add_table_column(label=title)
                for name, zi in comp.composition.items():
                    with dpg.table_row():
                        dpg.add_text(name)
                        dpg.add_text(f"{zi:.5f}")
                        for key, _title in _COMPOSITION_COLUMNS:
                            dpg.add_text(self._fmt(data.get(key, {}).get(name)))
            dpg.add_button(label="Close", callback=lambda: dpg.delete_item(win))

    def _on_delete_model_confirm(self, sender, app_data, user_data) -> None:
        mid = user_data
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        with dpg.window(label="Delete model", modal=True, no_resize=True,
                        width=420, height=140,
                        pos=(max(0, w // 2 - 210), max(0, h // 2 - 70))) as win:
            self._track_modal(win)
            dpg.add_text(f"Delete model '{mid}' from the database?")
            dpg.add_text("This cannot be undone.")
            dpg.add_spacer(height=8)
            with dpg.group(horizontal=True):
                dpg.add_button(label="Delete", width=120, user_data=(mid, win),
                               callback=self._on_delete_model_do)
                dpg.add_button(label="Cancel", width=120,
                               callback=lambda: dpg.delete_item(win))

    def _on_delete_model_do(self, sender, app_data, user_data) -> None:
        mid, win = user_data
        dpg.delete_item(win)
        try:
            ok = proj_svc.delete_model(self._state.db_path, mid)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Не удалось удалить модель %s", mid)
            self._set_status(f"Delete failed: {exc}")
            return
        if not ok:
            self._set_status(f"Model '{mid}' not found.")
            return
        # почистить след модели в сессии/выборе
        if self._session.workspaces:
            self._session.workspaces.pop(mid, None)
        if self._session.active_model_id == mid:
            self._session.active_model_id = None
        if self._selected_project == mid:
            self._selected_project = None
        self._state.refresh_model_list()  # notify -> перерисовка Projects
        self._set_status(f"Model '{mid}' deleted.")

    def _on_back_to_projects(self, sender=None, app_data=None, user_data=None) -> None:
        if self._active_job is not None:
            self._set_status("Calculation in progress - wait or cancel first.")
            return
        self._state.show_projects()

    def _on_new_fluid_manual(self, sender=None, app_data=None, user_data=None) -> None:
        self._new_fluid_form.clear_prefill()
        self._open_new_fluid_modal()

    def _open_new_fluid_modal(self) -> None:
        """Открывает форму нового флюида в модальном окне (не на весь экран)."""
        if self._new_fluid_win and dpg.does_item_exist(self._new_fluid_win):
            dpg.delete_item(self._new_fluid_win)
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        width, height = 680, max(420, h - 120)
        with dpg.window(label="New fluid", modal=True, no_collapse=True,
                        width=width, height=height,
                        pos=(max(0, w // 2 - width // 2), 60)) as win:
            self._new_fluid_win = self._track_modal(win)
            content = dpg.add_group()
        self._new_fluid_form.render(content)

    def _close_new_fluid_modal(self) -> None:
        if self._new_fluid_win and dpg.does_item_exist(self._new_fluid_win):
            dpg.delete_item(self._new_fluid_win)
        self._new_fluid_win = None

    def _on_fluid_created(self, model_id: str) -> None:
        """Колбэк формы: новая модель сохранена — закрыть модалку, открыть в workspace."""
        self._close_new_fluid_modal()
        self._state.refresh_model_list()
        self._restored_models.add(model_id)  # свежая модель, восстанавливать нечего
        self._expanded_models.add(model_id)
        self._state.enter_model(model_id)

    # --- импорт Excel -------------------------------------------------------

    def _on_import_excel(self, sender=None, app_data=None, user_data=None) -> None:
        """Открывает файловый диалог (создаётся по требованию, авто-id)."""
        dialog = dpg.add_file_dialog(
            label="Select Excel file with composition",
            directory_selector=False, show=True, modal=True,
            width=720, height=440,
            callback=self._on_excel_file_chosen,
            cancel_callback=lambda s, a: dpg.delete_item(s))
        self._track_modal(dialog)
        dpg.add_file_extension(".xlsx", parent=dialog)
        dpg.add_file_extension(".xls", parent=dialog)
        dpg.add_file_extension(".*", parent=dialog)

    def _on_excel_file_chosen(self, sender, app_data) -> None:
        if sender and dpg.does_item_exist(sender):
            dpg.delete_item(sender)
        path = (app_data or {}).get("file_path_name") or ""
        if not path:
            return
        self._excel_path = path
        try:
            sheets = proj_svc.excel_sheet_names(path)
        except Exception as exc:  # noqa: BLE001 — не Excel/битый файл
            logger.exception("Не удалось прочитать Excel %s", path)
            self._set_status(f"Cannot read Excel file: {exc}")
            return
        if not sheets:
            self._set_status("Excel file has no sheets.")
            return

        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        with dpg.window(label="Import composition from Excel", modal=True,
                        width=760, height=560,
                        pos=(max(0, w // 2 - 380), max(0, h // 2 - 280))) as win:
            self._excel_win = self._track_modal(win)
            dpg.add_text(Path(path).name)
            with dpg.group(horizontal=True):
                self._excel_sheet_id = dpg.add_combo(
                    items=sheets, default_value=sheets[0], width=220,
                    label="Sheet", callback=lambda: self._render_excel_preview())
                self._excel_header_id = dpg.add_checkbox(
                    label="First row is a header", default_value=True,
                    callback=lambda: self._render_excel_preview())
            dpg.add_text("Expected: two columns - component name | mole fraction.")
            self._excel_preview_group = dpg.add_group()
            dpg.add_separator()
            with dpg.group(horizontal=True):
                dpg.add_button(label="Load as composition", width=180,
                               callback=self._on_excel_load)
                dpg.add_button(label="Cancel", width=120,
                               callback=lambda: dpg.delete_item(win))
        self._render_excel_preview()

    def _render_excel_preview(self, *args) -> None:
        grp = self._excel_preview_group
        if grp is None or not dpg.does_item_exist(grp):
            return
        dpg.delete_item(grp, children_only=True)
        sheet = dpg.get_value(self._excel_sheet_id)
        header = bool(dpg.get_value(self._excel_header_id))
        try:
            prev = proj_svc.excel_preview(self._excel_path, sheet, header)
        except Exception as exc:  # noqa: BLE001
            msg = dpg.add_text(f"Cannot preview this sheet: {exc}", parent=grp)
            dpg.bind_item_theme(msg, self._theme_stale())
            return
        dpg.add_text(f"Preview ({sheet})"
                     + (" - first rows shown" if prev["truncated"] else ""),
                     parent=grp)
        with dpg.table(parent=grp, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       scrollY=True, scrollX=True, height=300, freeze_rows=1):
            for c in prev["columns"]:
                dpg.add_table_column(label=c)
            for row in prev["rows"]:
                with dpg.table_row():
                    for v in row:
                        dpg.add_text(v)

    def _on_excel_load(self, sender=None, app_data=None, user_data=None) -> None:
        sheet = dpg.get_value(self._excel_sheet_id)
        header = bool(dpg.get_value(self._excel_header_id))
        try:
            res = proj_svc.import_excel(self._excel_path, header, sheet)
        except Exception as exc:  # noqa: BLE001 — нечисловые значения/типы
            logger.exception("Импорт Excel не удался")
            self._set_status(f"Excel import failed: {exc}")
            return
        if not res["recognized"]:
            self._set_status("Excel import: no recognized components found.")
            return
        if self._excel_win and dpg.does_item_exist(self._excel_win):
            dpg.delete_item(self._excel_win)
        self._new_fluid_form.set_prefill(
            res["recognized"], name=Path(self._excel_path).stem,
            unrecognized=res["unrecognized"])
        self._open_new_fluid_modal()
        n_bad = len(res["unrecognized"])
        self._set_status(
            f"Imported {len(res['recognized'])} components from Excel"
            + (f", {n_bad} unrecognized skipped." if n_bad else "."))

    # ==================================================================
    #  Дерево (левая панель)
    # ==================================================================

    def _render_tree(self) -> None:
        """Дерево показывает ТОЛЬКО активную модель (выбор — на Projects)."""
        if not dpg.does_item_exist(_MODEL_TREE):
            return
        dpg.delete_item(_MODEL_TREE, children_only=True)

        model = self._state.active_model
        if model is None:
            dpg.add_text("No model selected.", parent=_MODEL_TREE)
            return
        expanded = model.model_id in self._expanded_models
        arrow = "v " if expanded else "> "
        dpg.add_selectable(
            label=f"{arrow}{model.title}  [{model.n_components}c, {model.eos}]",
            parent=_MODEL_TREE, default_value=True,
            user_data=model.model_id, callback=self._on_model_row,
        )
        if expanded and model.loaded:
            self._render_model_children(model)

    def _render_model_children(self, model) -> None:
        variant = model.variants.get("base")
        if variant is None:
            return
        is_active_model = model.model_id == self._state.active_model_id
        active_nid = variant.active_node_id if is_active_model else None

        # Composition
        comp = variant.nodes.get("composition")
        stale = "  *" if (comp and comp.status is NodeStatus.STALE) else ""
        dpg.add_selectable(
            label=f"    Composition{stale}", parent=_MODEL_TREE,
            default_value=(active_nid == "composition"),
            user_data=(model.model_id, "composition"),
            callback=self._on_tree_open_node,
        )

        # Flash (категория с историей запусков)
        cat_key = f"{model.model_id}:flash"
        cat_exp = cat_key in self._expanded_cats
        runs = variant.flash_runs()
        dpg.add_selectable(
            label=f"  {'v' if cat_exp else '>'} Flash ({len(runs)})",
            parent=_MODEL_TREE, user_data=cat_key, callback=self._on_cat_toggle,
        )
        if cat_exp:
            for run in runs:
                mark = "[*] " if run.node_id in self._compare_selection else ""
                sel = dpg.add_selectable(
                    label="      " + mark + self._flash_tree_label(run),
                    parent=_MODEL_TREE, default_value=(active_nid == run.node_id),
                    user_data=(model.model_id, run.node_id),
                    callback=self._on_tree_open_node,
                )
                self._attach_flash_context_menu(sel, run.node_id)
            dpg.add_selectable(label="      + New flash", parent=_MODEL_TREE,
                               user_data=model.model_id, callback=self._on_new_flash)
            n_sel = sum(1 for m in self._compare_selection if m in variant.nodes)
            if n_sel >= 2:
                dpg.add_selectable(label=f"      = Compare ({n_sel})",
                                   parent=_MODEL_TREE, user_data=model.model_id,
                                   callback=self._on_open_compare)

        # Experiments — вложенность по типам: Experiments → DLE/CCE/Separator → запуски
        mid = model.model_id
        ecat_key = f"{mid}:exp"
        ecat_exp = ecat_key in self._expanded_cats
        exp_runs = variant.experiment_runs()
        dpg.add_selectable(
            label=f"  {'v' if ecat_exp else '>'} Experiments ({len(exp_runs)})",
            parent=_MODEL_TREE, user_data=ecat_key, callback=self._on_cat_toggle)
        if ecat_exp:
            for kind in ("cce", "dle", "separator"):
                lbl = exp_svc.EXPERIMENT_TYPES[kind]["label"]
                kruns = [r for r in exp_runs if r.params.get("kind") == kind]
                kkey = f"{mid}:exp:{kind}"
                kexp = kkey in self._expanded_cats
                dpg.add_selectable(
                    label=f"    {'v' if kexp else '>'} {lbl} ({len(kruns)})",
                    parent=_MODEL_TREE, user_data=kkey, callback=self._on_cat_toggle)
                if kexp:
                    for run in kruns:
                        sel = dpg.add_selectable(
                            label="        " + self._exp_run_leaf_label(run),
                            parent=_MODEL_TREE,
                            default_value=(active_nid == run.node_id),
                            user_data=(mid, run.node_id),
                            callback=self._on_tree_open_node)
                        self._attach_exp_context_menu(sel, run.node_id)
                    dpg.add_selectable(label=f"        + New {lbl}",
                                       parent=_MODEL_TREE, user_data=(mid, kind),
                                       callback=self._on_new_experiment)

        dpg.add_text("    Saturation (soon)", parent=_MODEL_TREE)

    # --- обработчики дерева ----------------------------------------------

    def _on_model_row(self, sender, app_data, user_data) -> None:
        # модель уже активна (единственная в дереве) — только разворот
        mid = user_data
        if mid in self._expanded_models:
            self._expanded_models.discard(mid)
        else:
            self._expanded_models.add(mid)
        self._render_tree()

    def _on_cat_toggle(self, sender, app_data, user_data) -> None:
        cat_key = user_data
        if cat_key in self._expanded_cats:
            self._expanded_cats.discard(cat_key)
        else:
            self._expanded_cats.add(cat_key)
        self._render_tree()

    def _on_tree_open_node(self, sender, app_data, user_data) -> None:
        mid, nid = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid)
        node = self._state.node_by_id(nid)
        # Ctrl+клик по флэшу — переключить участие в сравнении
        if self._ctrl_down() and node is not None and node.kind is NodeKind.FLASH:
            self._toggle_compare(nid)
            return
        self._state.open_node(nid)

    def _toggle_compare(self, node_id: str) -> None:
        if node_id in self._compare_selection:
            self._compare_selection.discard(node_id)
        else:
            self._compare_selection.add(node_id)
        self._render_tree()
        self._set_status(f"Compare selection: {len(self._compare_selection)} run(s).")

    def _on_toggle_compare(self, sender, app_data, user_data) -> None:
        self._toggle_compare(user_data)

    def _on_open_compare(self, sender, app_data, user_data) -> None:
        mid = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid)
        members = [m for m in self._compare_selection
                   if m in self._state.active_variant.nodes]
        self._state.open_compare(members)

    def _on_new_flash(self, sender, app_data, user_data) -> None:
        mid = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid)
        self._state.new_flash_run()
        self._set_status("New flash tab opened - set P/T and Run.")

    def _attach_flash_context_menu(self, item_id, node_id: str) -> None:
        """Правый клик по листу флэша: переименовать / дублировать / удалить."""
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text("Flash run")
            dpg.add_separator()
            dpg.add_input_text(hint="rename + Enter", width=200, on_enter=True,
                               user_data=node_id, callback=self._on_flash_rename)
            dpg.add_button(label="Add / remove from compare", width=200,
                           user_data=node_id, callback=self._on_toggle_compare)
            dpg.add_button(label="Duplicate", width=200, user_data=node_id,
                           callback=self._on_flash_duplicate)
            dpg.add_button(label="Delete", width=200, user_data=node_id,
                           callback=self._on_flash_delete)

    def _on_flash_rename(self, sender, app_data, user_data) -> None:
        self._state.rename_node(user_data, app_data)
        self._set_status(f"Flash renamed to '{app_data}'." if app_data.strip()
                         else "Flash label cleared.")

    def _on_flash_duplicate(self, sender, app_data, user_data) -> None:
        self._state.duplicate_flash(user_data)
        self._set_status("Flash run duplicated.")

    def _on_flash_delete(self, sender, app_data, user_data) -> None:
        self._state.delete_node(user_data)
        self._set_status("Flash run deleted.")

    # --- эксперименты: дерево ---------------------------------------------

    def _exp_status_suffix(self, node) -> str:
        if node.status is NodeStatus.RUNNING:
            return " - running"
        if node.result is not None:
            return "" if node.status is NodeStatus.FRESH else " (stale)"
        return " - (not run)"

    def _exp_tree_label(self, node) -> str:
        """Подпись эксперимента с типом (для хлебных крошек)."""
        kind = node.params.get("kind", "")
        base = exp_svc.EXPERIMENT_TYPES.get(kind, {}).get("label", kind.upper())
        core = base + self._exp_status_suffix(node)
        label = node.params.get("label")
        return f"{label}  ({core})" if label else core

    def _exp_run_leaf_label(self, node) -> str:
        """Подпись листа под своей категорией — без повтора типа, с параметрами."""
        p = node.params
        seq = node.node_id.split("_")[-1]
        desc = f"T={self._g(p.get('T_c'))}C"
        if p.get("kind") in ("dle", "separator"):
            desc += f", Pres={self._g(p.get('P_res'))}"
        core = f"#{seq}  {desc}{self._exp_status_suffix(node)}"
        label = p.get("label")
        return f"{label}  ({core})" if label else core

    def _on_new_experiment(self, sender, app_data, user_data) -> None:
        mid, kind = user_data
        if mid != self._state.active_model_id:
            self._state.set_active_model(mid)
        composition = self._state.active_composition
        if composition is None:
            return
        pressures = exp_svc.default_pressures(composition)
        t_c = round(composition.T - 273.15, 2)
        defaults = {
            "pressures": pressures,
            "T_c": t_c,
            "P_res": max(pressures),
            "stage_temps_c": [15.0] * len(pressures),
        }
        self._state.new_experiment(kind, defaults)
        self._expanded_cats.add(f"{mid}:exp")
        self._expanded_cats.add(f"{mid}:exp:{kind}")
        self._set_status(f"New {kind.upper()} tab - set stages and Run.")

    def _attach_exp_context_menu(self, item_id, node_id: str) -> None:
        with dpg.popup(item_id, mousebutton=dpg.mvMouseButton_Right):
            dpg.add_text("Experiment")
            dpg.add_separator()
            dpg.add_input_text(hint="rename + Enter", width=200, on_enter=True,
                               user_data=node_id, callback=self._on_flash_rename)
            dpg.add_button(label="Duplicate", width=200, user_data=node_id,
                           callback=self._on_exp_duplicate)
            dpg.add_button(label="Delete", width=200, user_data=node_id,
                           callback=self._on_flash_delete)

    def _on_exp_duplicate(self, sender, app_data, user_data) -> None:
        node = self._state.node_by_id(user_data)
        if node is not None:
            defaults = {k: v for k, v in node.params.items()
                        if k not in ("kind", "label")}
            self._state.new_experiment(node.params["kind"], defaults)
            self._set_status("Experiment duplicated.")

    # ==================================================================
    #  Рабочая область (вкладки)
    # ==================================================================

    def _render_workspace(self) -> None:
        if not dpg.does_item_exist(_WORKSPACE):
            return
        dpg.delete_item(_WORKSPACE, children_only=True)
        self._flash_input_ids = {}
        self._exp_input_ids = {}
        self._exp_chart_holder = {}

        model = self._state.active_model
        variant = self._state.active_variant
        if model is None or variant is None:
            dpg.add_text("Select a model in the tree on the left.", parent=_WORKSPACE)
            return

        # хлебные крошки
        node = self._state.active_node
        crumb = model.title + (f"   >   {self._node_crumb(node)}" if node else "")
        dpg.add_text(crumb, parent=_WORKSPACE)
        dpg.add_separator(parent=_WORKSPACE)

        if not variant.open_node_ids:
            dpg.add_text("Open a node from the tree on the left.", parent=_WORKSPACE)
            return

        self._tab_ids = {}
        with dpg.tab_bar(parent=_WORKSPACE, reorderable=True,
                         callback=self._on_tab_changed) as tabbar:
            self._tabbar_id = tabbar
            for nid in variant.open_node_ids:
                n = variant.nodes.get(nid)
                if n is None:
                    continue
                with dpg.tab(label=self._tab_label(n)) as tab_id:
                    self._tab_ids[nid] = tab_id
                    page = dpg.add_group()
                    with dpg.group(horizontal=True, parent=page):
                        dpg.add_button(label="x Close", user_data=nid,
                                       callback=self._on_close_tab)
                    dpg.add_separator(parent=page)
                    if n.kind is NodeKind.COMPOSITION:
                        self._render_composition_tab(page, n)
                    elif n.kind is NodeKind.FLASH:
                        self._render_flash_tab(page, n)
                    elif n.kind is NodeKind.EXPERIMENT:
                        self._render_experiment_tab(page, n)
                    elif n.kind is NodeKind.COMPARE:
                        self._render_compare_tab(page, n)

        # синхронизировать активную вкладку
        active = variant.active_node_id
        if active in self._tab_ids:
            dpg.set_value(self._tabbar_id, self._tab_ids[active])

    def _on_tab_changed(self, sender, app_data, user_data=None) -> None:
        # app_data — id выбранной вкладки; обратное сопоставление к node_id
        nid = next((k for k, v in self._tab_ids.items() if v == app_data), None)
        if nid is not None:
            self._state.focus_node(nid)
            self._render_tree()  # обновить подсветку без пересборки вкладок

    def _on_close_tab(self, sender, app_data, user_data) -> None:
        self._state.close_node(user_data)

    def _tab_label(self, node) -> str:
        if node.kind is NodeKind.COMPOSITION:
            return "Composition"
        if node.kind is NodeKind.COMPARE:
            return f"Compare ({len(node.params.get('members', []))})"
        if node.kind is NodeKind.EXPERIMENT:
            label = node.params.get("label")
            base = exp_svc.EXPERIMENT_TYPES.get(node.params.get("kind"), {}).get(
                "label", "Exp")
            return label or f"{base} {node.node_id.split('_')[-1]}"
        if node.kind is NodeKind.FLASH:
            label = node.params.get("label")
            if label:
                return label
            p, t = node.params.get("P"), node.params.get("T")
            return f"Flash {self._g(p)}/{self._g(t)}"
        return node.title

    def _node_crumb(self, node) -> str:
        if node.kind is NodeKind.COMPOSITION:
            return "Composition"
        if node.kind is NodeKind.COMPARE:
            return f"Compare ({len(node.params.get('members', []))} runs)"
        if node.kind is NodeKind.EXPERIMENT:
            return "Experiments  >  " + self._exp_tree_label(node)
        if node.kind is NodeKind.FLASH:
            return "Flash  >  " + self._flash_tree_label(node)
        return node.title

    def _flash_tree_label(self, node) -> str:
        p, t = node.params.get("P"), node.params.get("T")
        if node.status is NodeStatus.RUNNING:
            core = f"{self._g(p)}/{self._g(t)} - running"
        elif node.result is not None:
            phase = "2-phase" if node.result.is_two_phase else "1-phase"
            suffix = "" if node.status is NodeStatus.FRESH else " (stale)"
            core = f"{self._g(p)} bar / {self._g(t)} C - {phase}{suffix}"
        else:
            core = f"{self._g(p)} bar / {self._g(t)} C - (not run)"
        label = node.params.get("label")
        return f"{label}  ({core})" if label else core

    # ==================================================================
    #  Вкладка Composition
    # ==================================================================

    def _render_composition_tab(self, parent, node) -> None:
        composition = self._state.active_composition
        if composition is None:
            dpg.add_text("No composition.", parent=parent)
            return

        status_line = dpg.add_text(f"Node status: {node.status.name}", parent=parent)
        if node.status is NodeStatus.STALE:
            dpg.bind_item_theme(status_line, self._theme_stale())
        if node.error:
            err = dpg.add_text(f"Error: {node.error}", parent=parent)
            dpg.bind_item_theme(err, self._theme_stale())
        dpg.add_separator(parent=parent)

        with dpg.tab_bar(parent=parent):
            with dpg.tab(label="Properties & correlations"):
                p = dpg.add_group()
                self._render_controls(composition, node, p)
                dpg.add_separator(parent=p)
                self._render_composition_table(composition, p)
            with dpg.tab(label="BIP matrix"):
                p = dpg.add_group()
                self._render_bip_matrix(composition, p)

    def _render_controls(self, composition, node, parent) -> None:
        params = node.params if node else {}
        eos_value = params.get("eos", composition.eos_name.value)
        correlations = params.get("correlations", comp_svc.default_correlations())

        with dpg.group(horizontal=True, parent=parent):
            dpg.add_combo(items=comp_svc.EOS_OPTIONS, default_value=eos_value,
                          width=140, label="EOS", callback=self._on_eos_changed)
            dpg.add_button(label="Recalculate properties",
                           callback=lambda: self._state.recalculate_composition())

        if comp_svc.has_c7_plus(composition):
            with dpg.tree_node(label="C7+ correlations", parent=parent,
                               default_open=True):
                for prop, prop_label in comp_svc.CORRELATION_LABELS:
                    dpg.add_combo(items=comp_svc.CORRELATION_OPTIONS[prop],
                                  default_value=correlations.get(prop, ""), width=220,
                                  label=prop_label, user_data=prop,
                                  callback=self._on_correlation_changed)
        else:
            dpg.add_text("No C7+ components - correlations not used.", parent=parent)

    def _render_composition_table(self, composition, parent) -> None:
        comp = composition.composition
        data = composition.composition_data

        with dpg.group(horizontal=True, parent=parent):
            dpg.add_text(f"Sum zi = {comp_svc.sum_zi(composition):.6f}")
            dpg.add_button(label="Normalize zi",
                           callback=lambda: self._state.normalize_composition())

        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, height=-1, freeze_rows=1):
            dpg.add_table_column(label="Component")
            dpg.add_table_column(label="zi, frac.")
            for _key, title in _COMPOSITION_COLUMNS:
                dpg.add_table_column(label=title)
            for name in comp:
                with dpg.table_row():
                    dpg.add_text(name)
                    dpg.add_input_float(default_value=float(comp[name]), width=90,
                                        step=0, format="%.5f", on_enter=True,
                                        user_data=name, callback=self._on_zi_edited)
                    for key, _title in _COMPOSITION_COLUMNS:
                        value = data.get(key, {}).get(name)
                        if value is None:
                            dpg.add_text("-")
                        else:
                            dpg.add_input_float(default_value=float(value), width=90,
                                                step=0, format="%.5g", on_enter=True,
                                                user_data=(name, key),
                                                callback=self._on_property_edited)

    def _render_bip_matrix(self, composition, parent) -> None:
        names = list(composition.composition.keys())
        data = composition.composition_data.get("bip", {})
        self._bip_ids = {}
        dpg.add_text("Symmetric matrix - editing cell (i, j) mirrors to (j, i). "
                     "Diagonal fixed at 0. Press Enter to apply.", parent=parent)
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       scrollX=True, scrollY=True, height=-1, freeze_rows=1,
                       freeze_columns=1, policy=dpg.mvTable_SizingFixedFit):
            dpg.add_table_column(label="")
            for name_j in names:
                dpg.add_table_column(label=name_j)
            for i, name_i in enumerate(names):
                with dpg.table_row():
                    dpg.add_text(name_i)
                    for j, name_j in enumerate(names):
                        if i == j:
                            dpg.add_input_float(default_value=0.0, width=_BIP_CELL_W,
                                                step=0, format="%.3f", readonly=True,
                                                enabled=False)
                        else:
                            value = float(data.get(name_i, {}).get(name_j, 0.0))
                            iid = dpg.add_input_float(
                                default_value=value, width=_BIP_CELL_W, step=0,
                                format="%.3f", on_enter=True, user_data=(i, j),
                                callback=self._on_bip_cell_edited)
                            self._bip_ids[(i, j)] = iid

    # --- обработчики Composition ----------------------------------------

    def _on_eos_changed(self, sender, app_data, user_data) -> None:
        self._state.set_composition_eos(app_data)
        self._set_status(f"EOS set to {app_data}.")

    def _on_correlation_changed(self, sender, app_data, user_data) -> None:
        self._state.set_correlation(user_data, app_data)
        self._set_status(f"Correlation {user_data} = '{app_data}' (press Recalculate).")

    def _on_zi_edited(self, sender, app_data, user_data) -> None:
        self._state.edit_zi(user_data, app_data)
        self._set_status(f"zi[{user_data}] = {app_data:.5f} (not normalized).")

    def _on_property_edited(self, sender, app_data, user_data) -> None:
        name, key = user_data
        self._state.edit_component_property(name, key, app_data)
        self._set_status(f"{key}[{name}] = {app_data:.5g}.")

    def _on_bip_cell_edited(self, sender, app_data, user_data) -> None:
        i, j = user_data
        composition = self._state.active_composition
        if composition is None:
            return
        names = list(composition.composition.keys())
        self._state.edit_bip(names[i], names[j], app_data)
        mirror = self._bip_ids.get((j, i))
        if mirror is not None and dpg.does_item_exist(mirror):
            dpg.set_value(mirror, app_data)
        self._set_status(f"BIP[{names[i]},{names[j]}] = {app_data:.5f}.")

    # ==================================================================
    #  Вкладка Flash
    # ==================================================================

    def _render_flash_tab(self, parent, node) -> None:
        nid = node.node_id
        params = node.params
        running = node.status is NodeStatus.RUNNING

        with dpg.group(horizontal=True, parent=parent):
            pid = dpg.add_input_float(label="P, bar", width=150, step=0,
                                      default_value=float(params.get("P", 100.0)),
                                      user_data=nid, callback=self._on_flash_param)
            tid = dpg.add_input_float(label="T, C", width=150, step=0,
                                      default_value=float(params.get("T", 20.0)),
                                      user_data=nid, callback=self._on_flash_param)
            self._flash_input_ids[nid] = (pid, tid)
            if running:
                dpg.add_loading_indicator(style=1, radius=2.0)
                dpg.add_button(label="Cancel", user_data=nid,
                               callback=self._on_flash_cancel)
            else:
                dpg.add_button(label="Run flash", user_data=nid,
                               callback=self._on_flash_run)

        if running:
            dpg.add_text("Running flash...", parent=parent)
        if node.status is NodeStatus.STALE and node.error:
            err = dpg.add_text(f"Flash error: {node.error}", parent=parent)
            dpg.bind_item_theme(err, self._theme_stale())
        elif node.status is NodeStatus.STALE and node.result is not None:
            st = dpg.add_text("Result is stale (composition changed) - re-run.",
                              parent=parent)
            dpg.bind_item_theme(st, self._theme_stale())

        if node.result is not None:
            self._render_flash_result(node.result, parent)

    def _render_flash_result(self, result, parent) -> None:
        dpg.add_separator(parent=parent)
        kind = "Two-phase" if result.is_two_phase else "Single-phase"
        dpg.add_text(f"{kind}    P = {result.pressure:.3f} bar    "
                     f"T = {result.temperature:.2f} K", parent=parent)
        # незакрываемые под-вкладки результата: Main + Composition
        with dpg.tab_bar(parent=parent):
            with dpg.tab(label="Main") as t_main:
                self._render_flash_main(result, t_main)
            with dpg.tab(label="Composition & K") as t_comp:
                self._render_flash_composition(result, t_comp)

    def _render_flash_main(self, result, parent) -> None:
        """Основные свойства фаз (доли, M, объём, плотность, Z, вязкость)."""
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, height=-1, freeze_rows=1):
            dpg.add_table_column(label="Property")
            dpg.add_table_column(label="Vapor")
            dpg.add_table_column(label="Liquid")
            with dpg.table_row():
                dpg.add_text("Phase mole fraction")
                dpg.add_text(self._fmt(result.vapor.mole_fraction))
                dpg.add_text(self._fmt(result.liquid.mole_fraction))
            for key, label in flash_service.PHASE_PROPERTY_ROWS:
                with dpg.table_row():
                    dpg.add_text(label)
                    dpg.add_text(self._fmt(result.vapor.properties.get(key)))
                    dpg.add_text(self._fmt(result.liquid.properties.get(key)))

    def _render_flash_composition(self, result, parent) -> None:
        """Состав каждой фазы (yi/xi) и константы равновесия K = yi/xi."""
        yi = result.vapor.composition if isinstance(result.vapor.composition, dict) else None
        xi = result.liquid.composition if isinstance(result.liquid.composition, dict) else None
        if not yi and not xi:
            dpg.add_text("Phase compositions are not available for this result.",
                         parent=parent)
            return

        names = list((xi or yi).keys())
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, height=-1, freeze_rows=1):
            dpg.add_table_column(label="Component")
            dpg.add_table_column(label="Vapor yi")
            dpg.add_table_column(label="Liquid xi")
            dpg.add_table_column(label="K = yi/xi")
            for name in names:
                y = yi.get(name) if yi else None
                x = xi.get(name) if xi else None
                k = (y / x) if (y is not None and x not in (None, 0.0)) else None
                with dpg.table_row():
                    dpg.add_text(name)
                    dpg.add_text(self._fmt(y))
                    dpg.add_text(self._fmt(x))
                    dpg.add_text(self._fmt(k))

    # ==================================================================
    #  Вкладка Compare (таблица сравнения + плитка панелей)
    # ==================================================================

    def _render_compare_tab(self, parent, node) -> None:
        members = [self._state.node_by_id(m) for m in node.params.get("members", [])]
        members = [m for m in members if m is not None and m.result is not None]
        if len(members) < 2:
            dpg.add_text("Select at least 2 computed flash runs to compare "
                         "(Ctrl+click runs in the tree, or use the right-click menu).",
                         parent=parent)
            return
        headers = [self._compare_col_label(m) for m in members]
        with dpg.tab_bar(parent=parent):
            with dpg.tab(label="Vapor") as t:
                self._compare_phase_table(members, headers, "vapor", t)
            with dpg.tab(label="Liquid") as t:
                self._compare_phase_table(members, headers, "liquid", t)
            with dpg.tab(label="K-values") as t:
                self._compare_k_table(members, headers, t)
            with dpg.tab(label="Panels") as t:
                self._compare_panels(members, headers, t)

    def _compare_col_label(self, node) -> str:
        label = node.params.get("label")
        if label:
            return label
        return f"{self._g(node.params.get('P'))}/{self._g(node.params.get('T'))}"

    def _compare_phase_table(self, members, headers, phase, parent) -> None:
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True, height=-1,
                       freeze_rows=1, freeze_columns=1):
            dpg.add_table_column(label="Property")
            for h in headers:
                dpg.add_table_column(label=h)
            with dpg.table_row():
                dpg.add_text("Phase mole fraction")
                for m in members:
                    dpg.add_text(self._fmt(getattr(m.result, phase).mole_fraction))
            for key, label in flash_service.PHASE_PROPERTY_ROWS:
                with dpg.table_row():
                    dpg.add_text(label)
                    for m in members:
                        dpg.add_text(self._fmt(getattr(m.result, phase).properties.get(key)))

    def _compare_k_table(self, members, headers, parent) -> None:
        # набор компонентов — из первого участника с известным составом жидкости
        names: list = []
        for m in members:
            xi = m.result.liquid.composition
            if isinstance(xi, dict):
                names = list(xi.keys())
                break
        if not names:
            dpg.add_text("Phase compositions are not available for these results.",
                         parent=parent)
            return
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True, height=-1,
                       freeze_rows=1, freeze_columns=1):
            dpg.add_table_column(label="Component (K=yi/xi)")
            for h in headers:
                dpg.add_table_column(label=h)
            for name in names:
                with dpg.table_row():
                    dpg.add_text(name)
                    for m in members:
                        yi = m.result.vapor.composition
                        xi = m.result.liquid.composition
                        y = yi.get(name) if isinstance(yi, dict) else None
                        x = xi.get(name) if isinstance(xi, dict) else None
                        k = (y / x) if (y is not None and x not in (None, 0.0)) else None
                        dpg.add_text(self._fmt(k))

    def _compare_panels(self, members, headers, parent) -> None:
        with dpg.group(horizontal=True, parent=parent):
            for m, h in zip(members, headers):
                with dpg.child_window(width=380, height=-1, border=True) as cw:
                    dpg.add_text(h)
                    self._render_flash_result(m.result, cw)

    def _flash_pt(self, nid: str) -> tuple[float, float]:
        """Текущие значения полей P/T вкладки флэша по её id узла."""
        pid, tid = self._flash_input_ids[nid]
        return dpg.get_value(pid), dpg.get_value(tid)

    # --- фоновая задача флэша (поток + опрос по кадрам) ------------------

    def _on_flash_param(self, sender, app_data, user_data) -> None:
        nid = user_data
        p, t = self._flash_pt(nid)
        self._state.set_flash_params(nid, p, t)

    def _on_flash_run(self, sender, app_data, user_data) -> None:
        if self._active_job is not None:
            self._set_status("Another flash is already running.")
            return
        nid = user_data
        composition = self._state.active_composition
        node = self._state.node_by_id(nid)
        if composition is None or node is None:
            return
        p, t = self._flash_pt(nid)
        self._state.set_flash_params(nid, p, t)

        holder = {"result": None, "error": None, "done": threading.Event()}
        self._active_job = {"holder": holder, "cancelled": False, "node_id": nid}
        self._state.set_flash_running(nid)

        def worker() -> None:
            try:
                holder["result"] = flash_service.run_flash(composition, p, t)
            except Exception as exc:  # noqa: BLE001
                holder["error"] = str(exc)
                logger.exception("Флэш-расчёт не удался")
            finally:
                holder["done"].set()

        threading.Thread(target=worker, daemon=True).start()
        self._set_status(f"Flash running at P={p:.3f} bar, T={t:.2f} C...")
        self._arm_flash_poll()

    def _arm_flash_poll(self) -> None:
        dpg.set_frame_callback(dpg.get_frame_count() + 1, self._poll_flash_job)

    def _poll_flash_job(self, sender=None, app_data=None) -> None:
        job = self._active_job
        if job is None:
            return
        if not job["holder"]["done"].is_set():
            self._arm_flash_poll()
            return
        self._active_job = None
        nid = job["node_id"]
        if job["cancelled"]:
            self._state.reset_flash(nid)
            self._set_status("Flash cancelled.")
        elif job["holder"]["error"]:
            self._state.set_flash_error(nid, job["holder"]["error"])
            self._set_status("Flash failed.")
        else:
            self._state.set_flash_result(nid, job["holder"]["result"])
            self._set_status("Flash done.")

    def _on_flash_cancel(self, sender, app_data, user_data) -> None:
        if self._active_job is not None:
            self._active_job["cancelled"] = True
            self._set_status("Cancelling flash (result will be discarded)...")

    # ==================================================================
    #  Вкладка Experiment (CCE / DLE / Separator)
    # ==================================================================

    def _render_experiment_tab(self, parent, node) -> None:
        nid = node.node_id
        p = node.params
        kind = p.get("kind", "cce")
        meta = exp_svc.EXPERIMENT_TYPES.get(kind, {})
        running = node.status is NodeStatus.RUNNING

        dpg.add_text(f"{meta.get('title', kind.upper())}", parent=parent)
        ids: dict = {}
        ids["pressures"] = dpg.add_input_text(
            label="Pressure stages, bar (comma-separated)", width=460, parent=parent,
            default_value=", ".join(self._g(x) for x in p.get("pressures", [])))
        with dpg.group(horizontal=True, parent=parent):
            ids["T_c"] = dpg.add_input_float(label="T, C", width=140, step=0,
                                             default_value=float(p.get("T_c", 100.0)))
            if meta.get("needs_p_res"):
                ids["P_res"] = dpg.add_input_float(
                    label="P res, bar", width=140, step=0,
                    default_value=float(p.get("P_res", 400.0)))
        if meta.get("needs_stage_temps"):
            ids["stage_temps"] = dpg.add_input_text(
                label="Stage T, C (comma-separated, same count)", width=460,
                parent=parent,
                default_value=", ".join(self._g(x) for x in p.get("stage_temps_c", [])))
        self._exp_input_ids[nid] = ids

        with dpg.group(horizontal=True, parent=parent):
            if running:
                dpg.add_loading_indicator(style=1, radius=2.0)
                dpg.add_button(label="Cancel", user_data=nid,
                               callback=self._on_flash_cancel)
            else:
                dpg.add_button(label="Run experiment", user_data=nid,
                               callback=self._on_experiment_run)
        if running:
            dpg.add_text("Running experiment (this may take a few seconds)...",
                         parent=parent)
        if node.status is NodeStatus.STALE and node.error:
            err = dpg.add_text(f"Error: {node.error}", parent=parent)
            dpg.bind_item_theme(err, self._theme_stale())
        elif node.status is NodeStatus.STALE and node.result is not None:
            st = dpg.add_text("Result is stale (composition changed) - re-run.",
                              parent=parent)
            dpg.bind_item_theme(st, self._theme_stale())

        if node.result is not None:
            dpg.add_separator(parent=parent)
            with dpg.tab_bar(parent=parent):
                with dpg.tab(label="Main") as t:
                    self._render_exp_table(node.result, t,
                                           node.result.get("main_columns", []))
                with dpg.tab(label="All data") as t:
                    self._render_exp_table(node.result, t, node.result["columns"])
                with dpg.tab(label="Composition by stage") as t:
                    self._render_exp_composition(node.result, t)
                with dpg.tab(label="Chart") as t:
                    self._render_exp_chart(node.result, t, nid)

    def _render_exp_table(self, result, parent, cols) -> None:
        allcols = result["columns"]
        idxs = [(c, allcols.index(c)) for c in cols if c in allcols]
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True, height=-1,
                       freeze_rows=1, freeze_columns=1):
            for c, _ in idxs:
                dpg.add_table_column(label=c)
            for row in result["rows"]:
                with dpg.table_row():
                    for _, i in idxs:
                        dpg.add_text(self._fmt(row[i]))

    def _render_exp_composition(self, result, parent) -> None:
        stages = result.get("stages", [])
        if not stages:
            dpg.add_text("Per-stage compositions are not available.", parent=parent)
            return
        with dpg.tab_bar(parent=parent):
            for phase in ("liquid", "vapor"):
                with dpg.tab(label=phase.capitalize()) as t:
                    self._render_exp_composition_phase(stages, phase, t)

    def _render_exp_composition_phase(self, stages, phase, parent) -> None:
        names: list = []
        for st in stages:
            if isinstance(st.get(phase), dict):
                names = list(st[phase].keys())
                break
        if not names:
            dpg.add_text(f"No {phase} composition.", parent=parent)
            return
        with dpg.table(parent=parent, header_row=True, borders_innerH=True,
                       borders_outerH=True, borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, scrollX=True, height=-1,
                       freeze_rows=1, freeze_columns=1):
            dpg.add_table_column(label="Component")
            for st in stages:
                dpg.add_table_column(label=f"{self._g(st.get('pressure'))} bar")
            for name in names:
                with dpg.table_row():
                    dpg.add_text(name)
                    for st in stages:
                        d = st.get(phase)
                        dpg.add_text(self._fmt(d.get(name) if isinstance(d, dict)
                                               else None))

    def _render_exp_chart(self, result, parent, nid) -> None:
        node = self._state.node_by_id(nid)
        charts = list(result.get("charts", []))
        for c in (node.params.get("extra_charts", []) if node else []):
            if c not in charts:
                charts.append(c)
        with dpg.group(horizontal=True, parent=parent):
            dpg.add_combo(items=result.get("plot_all", []), width=240,
                          label="Add chart for column", user_data=nid,
                          callback=self._on_exp_add_chart)
        holder = dpg.add_group(parent=parent)
        self._exp_chart_holder[nid] = holder
        if not charts:
            dpg.add_text("No plottable columns.", parent=holder)
        for col in charts:
            self._add_one_chart(result, col, holder)

    @staticmethod
    def _exp_x_range(result) -> tuple[float, float] | None:
        """Диапазон оси X (давления) по данным эксперимента, с малым отступом."""
        cols = result.get("columns", [])
        x = result.get("x")
        if x not in cols:
            return None
        xi = cols.index(x)
        xv = [row[xi] for row in result["rows"] if row[xi] is not None]
        if not xv:
            return None
        lo, hi = min(xv), max(xv)
        margin = (hi - lo) * 0.02 if hi > lo else 1.0
        return lo - margin, hi + margin

    def _add_one_chart(self, result, col, parent) -> None:
        xs, ys = exp_svc.series_for_plot(result, col)
        if not xs:
            return
        xr = self._exp_x_range(result)
        with dpg.plot(label=f"{col} vs pressure", height=220, width=-1, parent=parent):
            dpg.add_plot_legend()
            xax = dpg.add_plot_axis(dpg.mvXAxis, label="Pressure, bar")
            yax = dpg.add_plot_axis(dpg.mvYAxis, label=col)
            dpg.add_line_series(xs, ys, label=col, parent=yax)
            if xr is not None:
                dpg.set_axis_limits(xax, xr[0], xr[1])

    def _on_exp_add_chart(self, sender, app_data, user_data) -> None:
        nid = user_data
        col = app_data
        node = self._state.node_by_id(nid)
        if not col or node is None or node.result is None:
            return
        extra = node.params.setdefault("extra_charts", [])
        if col not in extra and col not in node.result.get("charts", []):
            extra.append(col)
        holder = self._exp_chart_holder.get(nid)
        if holder is not None and dpg.does_item_exist(holder):
            self._add_one_chart(node.result, col, holder)
        self._set_status(f"Added chart: {col}")

    def _on_experiment_run(self, sender, app_data, user_data) -> None:
        if self._active_job is not None:
            self._set_status("Another calculation is already running.")
            return
        nid = user_data
        node = self._state.node_by_id(nid)
        composition = self._state.active_composition
        if node is None or composition is None:
            return
        ids = self._exp_input_ids.get(nid, {})
        try:
            pressures = self._parse_floats(dpg.get_value(ids["pressures"]))
            if len(pressures) < 2:
                raise ValueError("need at least 2 pressure stages")
            params = dict(node.params)
            params["pressures"] = pressures
            params["T_c"] = float(dpg.get_value(ids["T_c"]))
            if "P_res" in ids:
                params["P_res"] = float(dpg.get_value(ids["P_res"]))
            if "stage_temps" in ids:
                temps = self._parse_floats(dpg.get_value(ids["stage_temps"]))
                if len(temps) != len(pressures):
                    raise ValueError("stage T count must match pressure count")
                params["stage_temps_c"] = temps
        except (ValueError, KeyError) as exc:
            self._state.set_node_error(nid, f"Invalid input: {exc}")
            return
        node.params.update(params)
        kind = node.params["kind"]

        holder = {"result": None, "error": None, "done": threading.Event()}
        self._active_job = {"holder": holder, "cancelled": False, "node_id": nid}
        self._state.set_node_running(nid)

        def worker() -> None:
            try:
                holder["result"] = exp_svc.run_experiment(composition, kind, params)
            except Exception as exc:  # noqa: BLE001
                holder["error"] = str(exc)
                logger.exception("Эксперимент %s не удался", kind)
            finally:
                holder["done"].set()

        threading.Thread(target=worker, daemon=True).start()
        self._set_status(f"{kind.upper()} running...")
        self._arm_flash_poll()

    @staticmethod
    def _parse_floats(text: str) -> list[float]:
        """Парсит числа из строки (разделители — запятая/пробел/перенос строки)."""
        out: list[float] = []
        for chunk in text.replace("\n", ",").replace(" ", ",").split(","):
            chunk = chunk.strip()
            if chunk:
                out.append(float(chunk))
        return out

    # ==================================================================
    #  Вспомогательное
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
        (флэши/эксперименты/вкладки). Повторный вызов для той же модели —
        no-op (guard `_restored_models`). Модель должна быть активной.
        """
        if model_id in self._restored_models:
            return
        self._restored_models.add(model_id)
        ws = (self._session.workspaces or {}).get(model_id) or {}
        variant = self._state.active_variant
        if variant is None or not ws:
            return

        # воссоздать флэш-узлы в исходном порядке (id flash_1, flash_2, … совпадут)
        for f in ws.get("flashes", []):
            nid = self._state.new_flash_run(f.get("P"), f.get("T"))
            snap = f.get("result")
            if nid and snap:
                try:
                    self._state.set_flash_result(
                        nid, flash_service.restore_flash_result(snap))
                except Exception:  # noqa: BLE001
                    logger.warning("Не удалось восстановить результат флэша")

        # воссоздать эксперименты (id exp_1, exp_2, … совпадут с сохранёнными)
        for e in ws.get("experiments", []):
            params = e.get("params", {})
            nid = self._state.new_experiment(
                e.get("kind"), {k: v for k, v in params.items() if k != "kind"})
            if nid and e.get("result"):
                self._state.set_node_result(nid, e["result"])

        # восстановить открытые вкладки / активную
        open_tabs = ws.get("open_tabs")
        if open_tabs is not None:
            variant.open_node_ids = [n for n in open_tabs if n in variant.nodes]
        active_tab = ws.get("active_tab")
        if active_tab and active_tab in variant.nodes:
            variant.active_node_id = active_tab
        elif variant.open_node_ids:
            variant.active_node_id = variant.open_node_ids[-1]

        # развернуть модель и категорию флэшей в дереве
        self._expanded_models.add(model_id)
        self._expanded_cats.add(f"{model_id}:flash")

        self._state.clear_history()  # восстановление не откатываем

    @staticmethod
    def _workspace_snapshot(variant) -> dict:
        """Сериализуемый снимок workspace варианта для сессии v2."""
        flashes: list = []
        for run in variant.flash_runs():
            item = {"P": run.params.get("P"), "T": run.params.get("T"),
                    "result": None}
            if run.result is not None and run.status is NodeStatus.FRESH:
                try:
                    item["result"] = flash_service.snapshot_flash_result(run.result)
                except Exception:  # noqa: BLE001
                    logger.warning("Не удалось сериализовать результат флэша")
            flashes.append(item)
        experiments: list = []
        for run in variant.experiment_runs():
            experiments.append({
                "kind": run.params.get("kind"),
                "params": dict(run.params),
                "result": run.result if run.status is NodeStatus.FRESH else None,
            })
        return {
            "open_tabs": list(variant.open_node_ids),
            "active_tab": variant.active_node_id,
            "flashes": flashes,
            "experiments": experiments,
        }

    def _persist_session(self) -> None:
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
