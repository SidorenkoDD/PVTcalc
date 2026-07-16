"""
Форма создания нового флюида (экран New fluid).

Особенности:
- Рендерится ОДИН раз при входе на экран (не по notify AppState) — весь ввод
  живёт в виджетах, полная перерисовка убила бы фокус ввода.
- Все id виджетов захватываются в словари (никаких фиксированных строковых
  тегов внутри перерисовываемого содержимого — правило проекта из-за
  отложенного освобождения алиасов в DPG).
- Компоненты — из каталога БД (`project_service.component_catalog`), тремя
  группами; у C7+ редактируемые M/gamma/Tb (дефолты БД). Сумма zi считается
  вживую через `set_value` (без перерисовки), есть Normalize.
- Создание — `project_service.create_model` (валидация, защита базы, T_res
  сохраняется новым полем записи); успех/отмена — через колбэки app.py.
- Поддерживает предзаполнение из Excel-импорта (`set_prefill`, Этап 4).
"""

import logging
from typing import Callable, Optional

import dearpygui.dearpygui as dpg

from gui.services import composition_service as comp_svc
from gui.services import project_service as proj_svc

logger = logging.getLogger(__name__)

_ZI_W = 100      # ширина поля zi
_OVR_W = 90      # ширина полей M/gamma/Tb
_TXT_W = 260     # ширина текстовых полей шапки


class NewFluidForm:
    """Форма «New fluid»; владелец экрана — PVTcalcApp (см. колбэки)."""

    def __init__(self, get_db_path: Callable[[], str],
                 on_created: Callable[[str], None],
                 on_cancel: Callable[[], None],
                 set_status: Callable[[str], None]):
        self._get_db_path = get_db_path
        self._on_created = on_created
        self._on_cancel = on_cancel
        self._set_status = set_status
        # захваченные id виджетов (пересоздаются при каждом render)
        self._ids: dict = {}
        self._zi_ids: dict[str, int] = {}
        self._ovr_ids: dict[str, dict] = {}
        self._corr_ids: dict[str, int] = {}
        # предзаполнение из импорта (Этап 4)
        self._prefill_zi: Optional[dict] = None
        self._prefill_name: Optional[str] = None
        self._unrecognized: Optional[dict] = None

    # --- предзаполнение (Excel) -------------------------------------------

    def set_prefill(self, zi: dict, name: Optional[str] = None,
                    unrecognized: Optional[dict] = None) -> None:
        """Задаёт zi/имя/список нераспознанных для следующего `render`."""
        self._prefill_zi = dict(zi)
        self._prefill_name = name
        self._unrecognized = dict(unrecognized) if unrecognized else None

    def clear_prefill(self) -> None:
        self._prefill_zi = None
        self._prefill_name = None
        self._unrecognized = None

    # --- отрисовка ----------------------------------------------------------

    def render(self, parent) -> None:
        dpg.delete_item(parent, children_only=True)
        self._ids = {}
        self._zi_ids = {}
        self._ovr_ids = {}
        self._corr_ids = {}

        dpg.add_spacer(height=6, parent=parent)
        with dpg.group(horizontal=True, parent=parent):
            dpg.add_text("New fluid")
            dpg.add_button(label="< Back to Projects",
                           callback=lambda: self._on_cancel())
        dpg.add_separator(parent=parent)

        if self._unrecognized:
            skipped = ", ".join(str(k) for k in self._unrecognized)
            warn = dpg.add_text(
                f"Unrecognized components skipped on import: {skipped}",
                parent=parent)
            dpg.bind_item_theme(warn, self._warn_theme())

        # --- шапка: имя/поле/EOS/T ---
        self._ids["name"] = dpg.add_input_text(
            label="Model name", width=_TXT_W, parent=parent,
            default_value=self._prefill_name or "")
        self._ids["model_id"] = dpg.add_input_text(
            label="Model id (empty = auto)", width=_TXT_W, parent=parent)
        self._ids["field"] = dpg.add_input_text(
            label="Field", width=_TXT_W, parent=parent)
        with dpg.group(horizontal=True, parent=parent):
            self._ids["eos"] = dpg.add_combo(
                items=comp_svc.EOS_OPTIONS, default_value="PREOS",
                width=140, label="EOS")
            self._ids["t_res"] = dpg.add_input_float(
                label="T res, K   (373.15 K = 100 C)", width=140, step=0,
                default_value=373.15)

        # --- сумма zi + действия ---
        with dpg.group(horizontal=True, parent=parent):
            self._ids["sum"] = dpg.add_text("Sum zi = 0.000000")
            dpg.add_button(label="Normalize", callback=self._on_normalize)
            dpg.add_button(label="Clear all", callback=self._on_clear)

        # --- компоненты по группам ---
        catalog = proj_svc.component_catalog()
        prefill = self._prefill_zi or {}
        for group, title in proj_svc.GROUP_TITLES.items():
            items = [c for c in catalog if c["group"] == group]
            with dpg.collapsing_header(label=title, parent=parent,
                                       default_open=True):
                self._render_group_table(items, group, prefill)

        # --- корреляции C7+ ---
        with dpg.collapsing_header(label="C7+ correlations", parent=parent,
                                   default_open=False):
            for prop, prop_label in comp_svc.CORRELATION_LABELS:
                self._corr_ids[prop] = dpg.add_combo(
                    items=comp_svc.CORRELATION_OPTIONS[prop],
                    default_value=comp_svc.DEFAULT_C7_CORRELATIONS[prop],
                    width=220, label=prop_label)

        # --- ошибки + создание ---
        self._ids["errors"] = dpg.add_text("", parent=parent, wrap=900)
        dpg.bind_item_theme(self._ids["errors"], self._warn_theme())
        with dpg.group(horizontal=True, parent=parent):
            dpg.add_button(label="Create fluid", callback=self._on_create)
            dpg.add_button(label="Cancel", callback=lambda: self._on_cancel())
        dpg.add_spacer(height=10, parent=parent)

        self._update_sum()

    def _render_group_table(self, items, group, prefill) -> None:
        is_c7 = group == proj_svc.GROUP_C7_PLUS
        with dpg.table(header_row=True, borders_innerH=True, borders_outerH=True,
                       borders_innerV=True, borders_outerV=True,
                       policy=dpg.mvTable_SizingFixedFit):
            dpg.add_table_column(label="Component")
            dpg.add_table_column(label="zi, frac.")
            if is_c7:
                dpg.add_table_column(label="M, g/mol")
                dpg.add_table_column(label="gamma")
                dpg.add_table_column(label="Tb, K")
            for c in items:
                name = c["name"]
                with dpg.table_row():
                    dpg.add_text(name)
                    self._zi_ids[name] = dpg.add_input_float(
                        default_value=float(prefill.get(name, 0.0)),
                        width=_ZI_W, step=0, format="%.6f",
                        callback=self._on_zi_changed)
                    if is_c7:
                        ovr = {}
                        for prop in ("molar_mass", "gamma", "Tb"):
                            ovr[prop] = dpg.add_input_float(
                                default_value=float(c.get(prop) or 0.0),
                                width=_OVR_W, step=0, format="%.4g")
                        self._ovr_ids[name] = ovr

    # --- обработчики ----------------------------------------------------------

    def _collect_zi(self) -> dict:
        return {name: dpg.get_value(iid) for name, iid in self._zi_ids.items()}

    def _update_sum(self) -> None:
        total = sum(v for v in self._collect_zi().values() if v)
        if dpg.does_item_exist(self._ids.get("sum", 0)):
            dpg.set_value(self._ids["sum"], f"Sum zi = {total:.6f}")

    def _on_zi_changed(self, sender=None, app_data=None, user_data=None) -> None:
        self._update_sum()

    def _on_normalize(self, sender=None, app_data=None, user_data=None) -> None:
        zi = self._collect_zi()
        total = sum(v for v in zi.values() if v and v > 0)
        if total <= 0:
            self._set_status("Nothing to normalize - all fractions are zero.")
            return
        for name, iid in self._zi_ids.items():
            v = zi.get(name) or 0.0
            dpg.set_value(iid, v / total if v > 0 else 0.0)
        self._update_sum()
        self._set_status("Fractions normalized to 1.")

    def _on_clear(self, sender=None, app_data=None, user_data=None) -> None:
        for iid in self._zi_ids.values():
            dpg.set_value(iid, 0.0)
        self._update_sum()

    def _show_errors(self, errors: list[str]) -> None:
        if dpg.does_item_exist(self._ids.get("errors", 0)):
            dpg.set_value(self._ids["errors"], "\n".join(errors))

    def _on_create(self, sender=None, app_data=None, user_data=None) -> None:
        db_path = self._get_db_path()
        name = dpg.get_value(self._ids["name"]).strip()
        model_id = dpg.get_value(self._ids["model_id"]).strip()
        if not model_id and name:
            model_id = proj_svc.suggest_model_id(name, db_path)
        zi = {k: v for k, v in self._collect_zi().items() if v and v > 0}

        errors = proj_svc.validate_new_model(model_id, name, zi, db_path)
        if errors:
            self._show_errors(errors)
            self._set_status("Fix the form errors and try again.")
            return
        self._show_errors([])

        correlations = {prop: dpg.get_value(cid)
                        for prop, cid in self._corr_ids.items()}
        c7_overrides = {
            comp_name: {prop: dpg.get_value(iid) for prop, iid in ovr.items()}
            for comp_name, ovr in self._ovr_ids.items() if comp_name in zi
        }
        try:
            proj_svc.create_model(
                db_path, model_id, name,
                dpg.get_value(self._ids["field"]).strip() or None,
                dpg.get_value(self._ids["eos"]),
                float(dpg.get_value(self._ids["t_res"])),
                zi, correlations, c7_overrides)
        except Exception as exc:  # noqa: BLE001 — показать любую ошибку движка
            logger.exception("Создание флюида не удалось")
            self._show_errors([f"Engine error: {exc}"])
            self._set_status("Fluid creation failed.")
            return

        self.clear_prefill()
        self._set_status(f"Fluid '{model_id}' created.")
        self._on_created(model_id)

    def _warn_theme(self):
        if getattr(self, "_warn_theme_id", None) is None:
            with dpg.theme() as theme_id:
                with dpg.theme_component(dpg.mvText):
                    dpg.add_theme_color(dpg.mvThemeCol_Text, (200, 110, 0))
            self._warn_theme_id = theme_id
        return self._warn_theme_id
