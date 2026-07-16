"""
Главное окно GUI на DearPyGui (тонкий View над `gui.app_state.AppState`).

Фаза 0+1 (вертикальный срез):
- одно главное окно на весь вьюпорт с **закреплённой раскладкой**
  (слева прибитая панель дерева, справа рабочая область, снизу статус) —
  без плавающего докинга, чтобы панели не «висели»;
- левая панель — дерево `Проект → Модель`, загрузка списка из `models.json`;
- рабочая область — таблица состава открытой модели (корневой узел);
- при выходе сохраняется высокоуровневая сессия (последняя модель + размер
  окна).

Весь видимый пользователю текст — на английском (комментарии/докстринги
остаются на русском по конвенции проекта). Узлы экспериментов (флэш и т.д.),
реактивная инвалидация и схема-граф — следующие фазы (см. docs/GUI.md).
"""

import logging
import threading
import time

import dearpygui.dearpygui as dpg

from gui.app_state import AppState, NodeStatus
from gui.services import composition_service as comp_svc
from gui.services import flash_service
from gui.session import SessionState, save_session
from gui.view.theme import build_light_theme

logger = logging.getLogger(__name__)

# Теги контейнеров, чтобы не пересоздавать окна и обновлять их содержимое
_PRIMARY = "primary_window"
_TREE_PANEL = "tree_panel"
_WORKSPACE_PANEL = "workspace_panel"
_MODEL_TREE = "model_tree"
_WORKSPACE = "workspace_content"
_STATUS_BAR = "status_bar"
_COMP_WINDOW = "composition_window"   # отдельное окно редактора состава
_COMP_CONTENT = "composition_content"
_RESTORE_MODAL = "restore_session_modal"

_STATUS_H = 28   # высота строки статуса внизу главного окна
_TREE_W = 320    # ширина закреплённой левой панели
_BIP_CELL_W = 60  # ширина ячейки матрицы BIP
_DOUBLE_CLICK_S = 0.4  # порог двойного клика по модели в дереве

# Скалярные свойства компонент, показываемые в таблице состава
# (ключ в composition_data -> заголовок колонки).
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
        # для распознавания двойного клика по модели в дереве
        self._last_click_model: str | None = None
        self._last_click_time: float = 0.0
        # активная фоновая задача флэша (расчёт в отдельном потоке)
        self._active_job: dict | None = None
        state.subscribe(self._render)

    # --- жизненный цикл --------------------------------------------------

    def run(self) -> None:
        dpg.create_context()
        dpg.bind_theme(build_light_theme())  # светлая тема на всё приложение

        dpg.create_viewport(
            title="PVTcalc",
            width=self._session.window_width,
            height=self._session.window_height,
        )

        self._build_menu()
        self._build_layout()

        dpg.setup_dearpygui()
        dpg.show_viewport()
        dpg.set_primary_window(_PRIMARY, True)

        # Начальное наполнение
        self._state.refresh_model_list()
        self._maybe_prompt_restore()

        dpg.start_dearpygui()

        # Сохранение сессии на выходе
        self._persist_session()
        dpg.destroy_context()

    # --- построение статичного каркаса ----------------------------------

    def _build_menu(self) -> None:
        with dpg.viewport_menu_bar():
            with dpg.menu(label="Project"):
                dpg.add_menu_item(label="Refresh models",
                                  callback=lambda: self._state.refresh_model_list())

    def _build_layout(self) -> None:
        """
        Главное окно на весь вьюпорт с фиксированной раскладкой:
        [ tree | workspace ] сверху, строка статуса снизу.
        """
        with dpg.window(tag=_PRIMARY, no_title_bar=True, no_move=True,
                        no_resize=True, no_collapse=True):
            with dpg.group(horizontal=True):
                # Левая панель дерева — прибита к левому краю фиксированной шириной
                with dpg.child_window(tag=_TREE_PANEL, width=_TREE_W,
                                      height=-_STATUS_H, border=True):
                    dpg.add_text("Models (models.json)")
                    dpg.add_separator()
                    dpg.add_group(tag=_MODEL_TREE)
                # Рабочая область — занимает остаток по ширине и высоте
                with dpg.child_window(tag=_WORKSPACE_PANEL, width=-1,
                                      height=-_STATUS_H, border=True):
                    dpg.add_group(tag=_WORKSPACE)
            dpg.add_separator()
            dpg.add_text("Ready.", tag=_STATUS_BAR)

    # --- реактивная перерисовка -----------------------------------------

    def _render(self) -> None:
        """Полная перерисовка динамических панелей по текущему `AppState`."""
        self._render_tree()
        self._render_workspace()
        self._render_composition_content()

    def _render_tree(self) -> None:
        if not dpg.does_item_exist(_MODEL_TREE):
            return
        dpg.delete_item(_MODEL_TREE, children_only=True)
        for model in self._state.models.values():
            label = f"{model.title}  [{model.n_components} comp., {model.eos}]"
            is_active = model.model_id == self._state.active_model_id
            dpg.add_selectable(
                label=("> " if is_active else "  ") + label,
                parent=_MODEL_TREE,
                callback=self._on_model_selected,
                user_data=model.model_id,
            )

    def _render_workspace(self) -> None:
        """Центральная область — эксперименты активной модели (флэш и т.д.)."""
        if not dpg.does_item_exist(_WORKSPACE):
            return
        dpg.delete_item(_WORKSPACE, children_only=True)

        model = self._state.active_model
        if model is None:
            dpg.add_text("Select a model in the tree on the left "
                         "(double-click to edit its composition).", parent=_WORKSPACE)
            return

        dpg.add_text(f"Model: {model.title}    Field: {model.field_name or '-'}",
                     parent=_WORKSPACE)
        node = self._state.composition_node
        status = node.status.name if node else "-"
        dpg.add_text(f"Composition: {status}   "
                     "(double-click the model in the tree to edit it)",
                     parent=_WORKSPACE)
        dpg.add_separator(parent=_WORKSPACE)

        dpg.add_text("Experiments", parent=_WORKSPACE)
        self._render_flash_panel()

    # --- панель флэша (центральная область) -----------------------------

    def _render_flash_panel(self) -> None:
        node = self._state.flash_node
        if node is None:
            return
        params = node.params
        running = node.status is NodeStatus.RUNNING

        with dpg.group(horizontal=True, parent=_WORKSPACE):
            dpg.add_text("Flash:")
            dpg.add_input_float(tag="flash_p", label="P, bar", width=150, step=0,
                                default_value=float(params.get("P", 100.0)),
                                callback=self._on_flash_param)
            dpg.add_input_float(tag="flash_t", label="T, C", width=150, step=0,
                                default_value=float(params.get("T", 100.0)),
                                callback=self._on_flash_param)
            if running:
                dpg.add_loading_indicator(style=1, radius=2.0)
                dpg.add_button(label="Cancel", callback=self._on_flash_cancel)
            else:
                dpg.add_button(label="Run flash", callback=self._on_flash_run)

        if running:
            dpg.add_text("Running flash...", parent=_WORKSPACE)

        if node.status is NodeStatus.STALE and node.error:
            err = dpg.add_text(f"Flash error: {node.error}", parent=_WORKSPACE)
            dpg.bind_item_theme(err, self._theme_stale())
        elif node.status is NodeStatus.STALE and node.result is not None:
            stale = dpg.add_text("Result is stale (composition changed) - re-run.",
                                 parent=_WORKSPACE)
            dpg.bind_item_theme(stale, self._theme_stale())

        if node.result is not None:
            self._render_flash_result(node.result)

    def _render_flash_result(self, result) -> None:
        dpg.add_separator(parent=_WORKSPACE)
        kind = "Two-phase" if result.is_two_phase else "Single-phase"
        dpg.add_text(f"{kind}    P = {result.pressure:.3f} bar    "
                     f"T = {result.temperature:.2f} K", parent=_WORKSPACE)

        with dpg.table(parent=_WORKSPACE, header_row=True,
                       borders_innerH=True, borders_outerH=True,
                       borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, height=-1):
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

    @staticmethod
    def _fmt(value) -> str:
        if value is None:
            return "-"
        if isinstance(value, (int, float)):
            return f"{value:.5g}"
        return str(value)

    # --- фоновая задача флэша (поток + опрос по кадрам) ------------------

    def _on_flash_param(self, sender, app_data, user_data) -> None:
        # синхронизируем введённые P/T в параметры узла (без перерисовки)
        p = dpg.get_value("flash_p")
        t = dpg.get_value("flash_t")
        self._state.set_flash_params(p, t)

    def _on_flash_run(self, sender, app_data, user_data) -> None:
        if self._active_job is not None:
            return
        composition = self._state.active_composition
        node = self._state.flash_node
        if composition is None or node is None:
            return
        p = dpg.get_value("flash_p")
        t = dpg.get_value("flash_t")
        self._state.set_flash_params(p, t)

        holder = {"result": None, "error": None, "done": threading.Event()}
        self._active_job = {"holder": holder, "cancelled": False}
        self._state.set_flash_running()  # notify -> перерисовка с индикатором

        def worker() -> None:
            try:
                holder["result"] = flash_service.run_flash(composition, p, t)
            except Exception as exc:  # noqa: BLE001 — донести ошибку в UI
                holder["error"] = str(exc)
                logger.exception("Флэш-расчёт не удался")
            finally:
                holder["done"].set()

        threading.Thread(target=worker, daemon=True).start()
        self._set_status(f"Flash running at P={p:.3f} bar, T={t:.2f} C...")
        self._arm_flash_poll()

    def _arm_flash_poll(self) -> None:
        # DPG-safe: результат применяется в главном потоке в frame-callback
        dpg.set_frame_callback(dpg.get_frame_count() + 1, self._poll_flash_job)

    def _poll_flash_job(self, sender=None, app_data=None) -> None:
        job = self._active_job
        if job is None:
            return
        if not job["holder"]["done"].is_set():
            self._arm_flash_poll()
            return

        self._active_job = None
        if job["cancelled"]:
            self._state.clear_flash()
            self._set_status("Flash cancelled.")
            return
        if job["holder"]["error"]:
            self._state.set_flash_error(job["holder"]["error"])
            self._set_status("Flash failed.")
        else:
            self._state.set_flash_result(job["holder"]["result"])
            self._set_status("Flash done.")

    def _on_flash_cancel(self, sender, app_data, user_data) -> None:
        # Движок не прерываемый — отменяем на уровне UI: результат отбрасывается,
        # когда поток завершится (см. `_poll_flash_job`).
        if self._active_job is not None:
            self._active_job["cancelled"] = True
            self._set_status("Cancelling flash (result will be discarded)...")

    # --- отдельное окно редактора состава -------------------------------

    def _open_composition_window(self, model_id: str) -> None:
        """Открывает (или показывает поверх) отдельное окно редактора состава."""
        if self._state.active_model_id != model_id:
            self._state.open_model(model_id)  # активирует + notify (окна ещё нет)
        if not dpg.does_item_exist(_COMP_WINDOW):
            with dpg.window(tag=_COMP_WINDOW, label="Composition", width=940,
                            height=680, pos=(380, 60), show=True):
                dpg.add_group(tag=_COMP_CONTENT)
        else:
            dpg.configure_item(_COMP_WINDOW, show=True)
        dpg.focus_item(_COMP_WINDOW)
        self._render_composition_content()

    def _render_composition_content(self) -> None:
        """Перерисовывает содержимое окна редактора состава, если оно открыто."""
        if not dpg.does_item_exist(_COMP_CONTENT):
            return
        if not dpg.is_item_shown(_COMP_WINDOW):
            return
        dpg.delete_item(_COMP_CONTENT, children_only=True)

        composition = self._state.active_composition
        model = self._state.active_model
        node = self._state.composition_node
        if composition is None or model is None:
            dpg.add_text("No active model.", parent=_COMP_CONTENT)
            return

        dpg.configure_item(_COMP_WINDOW, label=f"Composition - {model.title}")

        status = node.status if node else NodeStatus.EMPTY
        status_line = dpg.add_text(f"Node status: {status.name}", parent=_COMP_CONTENT)
        if status is NodeStatus.STALE:
            dpg.bind_item_theme(status_line, self._theme_stale())
        if node and node.error:
            err = dpg.add_text(f"Error: {node.error}", parent=_COMP_CONTENT)
            dpg.bind_item_theme(err, self._theme_stale())
        dpg.add_separator(parent=_COMP_CONTENT)

        with dpg.tab_bar(parent=_COMP_CONTENT):
            with dpg.tab(label="Properties & correlations"):
                page = dpg.add_group()
                self._render_controls(composition, node, page)
                dpg.add_separator(parent=page)
                self._render_composition_table(composition, page)
            with dpg.tab(label="BIP matrix"):
                page = dpg.add_group()
                self._render_bip_matrix(composition, page)

    # --- панель 1: EOS + корреляции + таблица состава -------------------

    def _render_controls(self, composition, node, parent) -> None:
        params = node.params if node else {}
        eos_value = params.get("eos", composition.eos_name.value)
        correlations = params.get("correlations", comp_svc.default_correlations())

        with dpg.group(horizontal=True, parent=parent):
            dpg.add_combo(
                items=comp_svc.EOS_OPTIONS, default_value=eos_value, width=140,
                label="EOS", callback=self._on_eos_changed,
            )
            dpg.add_button(label="Recalculate properties",
                           callback=lambda: self._state.recalculate_composition())

        if comp_svc.has_c7_plus(composition):
            with dpg.tree_node(label="C7+ correlations", parent=parent,
                               default_open=True):
                for prop, prop_label in comp_svc.CORRELATION_LABELS:
                    dpg.add_combo(
                        items=comp_svc.CORRELATION_OPTIONS[prop],
                        default_value=correlations.get(prop, ""),
                        width=220, label=prop_label,
                        user_data=prop, callback=self._on_correlation_changed,
                    )
        else:
            dpg.add_text("No C7+ components - correlations not used.", parent=parent)

    def _render_composition_table(self, composition, parent) -> None:
        comp = composition.composition
        data = composition.composition_data
        zi_sum = comp_svc.sum_zi(composition)

        with dpg.group(horizontal=True, parent=parent):
            dpg.add_text(f"Sum zi = {zi_sum:.6f}")
            dpg.add_button(label="Normalize zi",
                           callback=lambda: self._state.normalize_composition())

        with dpg.table(parent=parent, header_row=True,
                       borders_innerH=True, borders_outerH=True,
                       borders_innerV=True, borders_outerV=True,
                       resizable=True, scrollY=True, height=-1):
            dpg.add_table_column(label="Component")
            dpg.add_table_column(label="zi, frac.")
            for _key, title in _COMPOSITION_COLUMNS:
                dpg.add_table_column(label=title)

            for name in comp:
                with dpg.table_row():
                    dpg.add_text(name)
                    dpg.add_input_float(
                        default_value=float(comp[name]), width=90, step=0,
                        format="%.5f", on_enter=True,
                        user_data=name, callback=self._on_zi_edited,
                    )
                    for key, _title in _COMPOSITION_COLUMNS:
                        value = data.get(key, {}).get(name)
                        if value is None:
                            dpg.add_text("-")
                        else:
                            dpg.add_input_float(
                                default_value=float(value), width=90, step=0,
                                format="%.5g", on_enter=True,
                                user_data=(name, key), callback=self._on_property_edited,
                            )

    # --- панель 2: матрица BIP (симметричная, редактируемая) ------------

    def _render_bip_matrix(self, composition, parent) -> None:
        names = list(composition.composition.keys())
        data = composition.composition_data.get("bip", {})

        dpg.add_text("Symmetric matrix — editing cell (i, j) mirrors to (j, i). "
                     "Diagonal is fixed at 0. Press Enter to apply.", parent=parent)

        with dpg.table(parent=parent, header_row=True,
                       borders_innerH=True, borders_outerH=True,
                       borders_innerV=True, borders_outerV=True,
                       scrollX=True, scrollY=True, height=-1,
                       freeze_rows=1, freeze_columns=1,
                       policy=dpg.mvTable_SizingFixedFit):
            dpg.add_table_column(label="")  # угол: имена строк
            for name_j in names:
                dpg.add_table_column(label=name_j)

            for i, name_i in enumerate(names):
                with dpg.table_row():
                    dpg.add_text(name_i)
                    for j, name_j in enumerate(names):
                        if i == j:
                            dpg.add_input_float(default_value=0.0, width=_BIP_CELL_W,
                                                step=0, format="%.3f",
                                                readonly=True, enabled=False)
                        else:
                            value = float(data.get(name_i, {}).get(name_j, 0.0))
                            dpg.add_input_float(
                                default_value=value, width=_BIP_CELL_W, step=0,
                                format="%.3f", on_enter=True,
                                tag=self._bip_tag(i, j),
                                user_data=(i, j), callback=self._on_bip_cell_edited,
                            )

    @staticmethod
    def _bip_tag(i: int, j: int) -> str:
        return f"bipm_{i}_{j}"

    # --- обработчики Фазы 2 ---------------------------------------------

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
        # зеркалим значение в симметричную ячейку без полной перерисовки
        mirror = self._bip_tag(j, i)
        if dpg.does_item_exist(mirror):
            dpg.set_value(mirror, app_data)
        self._set_status(f"BIP[{names[i]},{names[j]}] = {app_data:.5f}.")

    def _theme_stale(self):
        """Ленивая тема (тёмно-оранжевый текст) для пометки устаревшего/ошибочного."""
        if getattr(self, "_stale_theme_id", None) is None:
            with dpg.theme() as theme_id:
                with dpg.theme_component(dpg.mvText):
                    dpg.add_theme_color(dpg.mvThemeCol_Text, (200, 110, 0))
            self._stale_theme_id = theme_id
        return self._stale_theme_id

    # --- обработчики -----------------------------------------------------

    def _on_model_selected(self, sender, app_data, user_data) -> None:
        """Один клик — выбрать модель; двойной — открыть окно редактора состава."""
        model_id = user_data
        now = time.monotonic()
        is_double = (self._last_click_model == model_id
                     and (now - self._last_click_time) < _DOUBLE_CLICK_S)
        self._last_click_model = None if is_double else model_id
        self._last_click_time = now
        try:
            if is_double:
                self._open_composition_window(model_id)
                self._set_status(f"Composition editor opened for '{model_id}'.")
            else:
                self._state.open_model(model_id)
                self._set_status(f"Model '{model_id}' selected.")
        except Exception as exc:  # noqa: BLE001 — показать пользователю любую ошибку загрузки
            logger.exception("Не удалось открыть модель %s", model_id)
            self._set_status(f"Failed to load '{model_id}': {exc}")

    def _set_status(self, text: str) -> None:
        if dpg.does_item_exist(_STATUS_BAR):
            dpg.set_value(_STATUS_BAR, text)

    # --- сессия ----------------------------------------------------------

    def _maybe_prompt_restore(self) -> None:
        """
        Если есть восстановимая сессия — показать модальный вопрос
        «продолжить прошлую сессию?»; иначе тихо стартовать с чистого листа.
        """
        mid = self._session.active_model_id
        if not mid or mid not in self._state.models:
            return

        date_str = self._format_saved_at(self._session.saved_at)
        w, h = self._session.window_width, self._session.window_height
        with dpg.window(label="Restore session", modal=True, no_collapse=True,
                        no_resize=True, no_move=True, tag=_RESTORE_MODAL,
                        width=460, height=170,
                        pos=(max(0, w // 2 - 230), max(0, h // 2 - 85))):
            if date_str:
                dpg.add_text(f"A previous session was saved on {date_str}.")
            else:
                dpg.add_text("A previous session was found.")
            dpg.add_text(f"Model: {mid}")
            dpg.add_spacer(height=10)
            dpg.add_text("Continue where you left off?")
            dpg.add_spacer(height=8)
            with dpg.group(horizontal=True):
                dpg.add_button(label="Continue", width=130,
                               callback=self._on_restore_continue)
                dpg.add_button(label="Start fresh", width=130,
                               callback=self._on_restore_decline)

    def _on_restore_continue(self, sender=None, app_data=None) -> None:
        if dpg.does_item_exist(_RESTORE_MODAL):
            dpg.delete_item(_RESTORE_MODAL)
        self._restore_session()

    def _on_restore_decline(self, sender=None, app_data=None) -> None:
        if dpg.does_item_exist(_RESTORE_MODAL):
            dpg.delete_item(_RESTORE_MODAL)
        self._set_status("Started fresh (previous session kept on disk).")

    @staticmethod
    def _format_saved_at(iso: str | None) -> str:
        """ISO-строку `saved_at` -> человекочитаемое 'YYYY-MM-DD HH:MM' (или '')."""
        if not iso:
            return ""
        try:
            from datetime import datetime
            return datetime.fromisoformat(iso).strftime("%Y-%m-%d %H:%M")
        except (ValueError, TypeError):
            return ""

    def _restore_session(self) -> None:
        """Восстанавливает прошлую сессию: модель, параметры/результат флэша,
        открытое окно редактора состава."""
        mid = self._session.active_model_id
        if not mid or mid not in self._state.models:
            return
        try:
            self._state.open_model(mid)
        except Exception:  # noqa: BLE001
            logger.warning("Не удалось восстановить модель сессии %s", mid)
            return

        # параметры + результат флэша
        flash = self._session.flash
        node = self._state.flash_node
        if flash and node is not None:
            if flash.get("P") is not None:
                node.params["P"] = flash["P"]
            if flash.get("T") is not None:
                node.params["T"] = flash["T"]
            snap = flash.get("result")
            if snap:
                try:
                    node.result = flash_service.restore_flash_result(snap)
                    node.status = NodeStatus.FRESH
                except Exception:  # noqa: BLE001 — битый слепок не должен ронять старт
                    logger.warning("Не удалось восстановить результат флэша из сессии")

        # окно редактора состава
        if self._session.composition_window_open:
            self._open_composition_window(mid)

        self._state._notify()  # перерисовать с восстановленным состоянием
        self._set_status(f"Session restored: model '{mid}'.")

    def _persist_session(self) -> None:
        self._session.active_model_id = self._state.active_model_id
        try:
            self._session.window_width = dpg.get_viewport_width()
            self._session.window_height = dpg.get_viewport_height()
        except Exception:  # noqa: BLE001 — viewport уже уничтожен
            pass

        self._session.composition_window_open = (
            dpg.does_item_exist(_COMP_WINDOW) and dpg.is_item_shown(_COMP_WINDOW)
        )

        node = self._state.flash_node
        if node is not None:
            flash: dict = {"P": node.params.get("P"), "T": node.params.get("T"),
                           "result": None}
            if node.result is not None and node.status is NodeStatus.FRESH:
                try:
                    flash["result"] = flash_service.snapshot_flash_result(node.result)
                except Exception:  # noqa: BLE001 — не мешать сохранению остального
                    logger.warning("Не удалось сериализовать результат флэша")
            self._session.flash = flash

        save_session(self._session)
