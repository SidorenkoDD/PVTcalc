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
import time

import dearpygui.dearpygui as dpg

from gui.app_state import AppState, NodeStatus
from gui.services import composition_service as comp_svc
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
        self._restore_session_selection()

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
        dpg.add_text("Flash (P, T) - coming in the next phase.", parent=_WORKSPACE)

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

    def _restore_session_selection(self) -> None:
        mid = self._session.active_model_id
        if mid and mid in self._state.models:
            try:
                self._state.open_model(mid)
            except Exception:  # noqa: BLE001
                logger.warning("Не удалось восстановить модель сессии %s", mid)

    def _persist_session(self) -> None:
        self._session.active_model_id = self._state.active_model_id
        try:
            self._session.window_width = dpg.get_viewport_width()
            self._session.window_height = dpg.get_viewport_height()
        except Exception:  # noqa: BLE001 — viewport уже уничтожен
            pass
        save_session(self._session)
