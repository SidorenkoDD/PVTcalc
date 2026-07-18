"""Страница Projects и диалоги создания/импорта состава.

Mixin не владеет состоянием приложения: он использует публичный контракт
composition root `PVTcalcApp` и изолирует крупный DPG-срез представления.
"""

import logging
import time
from pathlib import Path

import dearpygui.dearpygui as dpg

from gui.app_state import StateChange, StateChangeKind
from gui.services import project_service as proj_svc
from gui.view.contracts import ContextBoundView
from gui.view.new_fluid_form import NewFluidForm

logger = logging.getLogger(__name__)

_PROJECTS_CONTENT = "projects_content"

# Колонки read-only предпросмотра состава модели на Projects.
_COMPOSITION_COLUMNS: list[tuple[str, str]] = [
    ("molar_mass", "M, g/mol"),
    ("critical_temperature", "Tc, K"),
    ("critical_pressure", "Pc, bar"),
    ("acentric_factor", "omega"),
    ("critical_volume", "Vc"),
    ("shift_parameter", "shift"),
    ("Kw", "Kw"),
]


class ProjectsViewMixin(ContextBoundView):
    """Рендер Projects и callbacks создания/импорта моделей."""

    _new_fluid_form: NewFluidForm
    _new_fluid_win: int | None
    _selected_project: str | None
    _proj_row_ids: dict[str, int]
    _last_proj_click: str | None
    _last_proj_click_time: float
    _expanded_models: set[str]
    _restored_models: set[str]
    _excel_path: str | None
    _excel_win: int | None
    _excel_sheet_id: int | None
    _excel_header_id: int | None
    _excel_preview_group: int | None
    _e300_path: str | None
    _e300_win: int | None
    _duplicate_win: int | None
    _duplicate_ids: dict[str, int]

    def _render_projects(self) -> None:
        if not dpg.does_item_exist(_PROJECTS_CONTENT):
            return
        dpg.delete_item(_PROJECTS_CONTENT, children_only=True)
        parent = _PROJECTS_CONTENT

        dpg.add_spacer(height=6, parent=parent)
        dpg.add_text("Projects", parent=parent)
        dpg.add_separator(parent=parent)

        if self._state.model_list_error:
            warning = dpg.add_text(
                "Models database could not be read. No data was changed.\n"
                f"{self._state.model_list_error}",
                parent=parent,
                wrap=920,
            )
            dpg.bind_item_theme(warning, self._theme_stale())
            with dpg.group(horizontal=True, parent=parent):
                dpg.add_text(
                    "Restore or correct models.json, then refresh the list.",
                )
                dpg.add_button(
                    label="Refresh models",
                    callback=lambda: self._state.refresh_model_list(),
                )
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
                        label=model.title + ("  *" if model.dirty else ""), span_columns=True,
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
            dpg.add_button(label="Import E300", callback=self._on_import_e300)

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
        if self._jobs.busy:
            self._set_status("Calculation in progress - wait or cancel first.")
            return
        try:
            self._state.enter_model(mid, notify=False)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Не удалось открыть модель %s", mid)
            self._set_status(f"Failed to open '{mid}': {exc}")
            return
        self._expanded_models.add(mid)
        self._restore_workspace(mid)
        self._state.notify(StateChange(StateChangeKind.NAVIGATION))
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
            dpg.add_menu_item(label="Duplicate model...", user_data=model_id,
                              callback=self._on_duplicate_model_confirm)
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

    def _on_duplicate_model_confirm(self, sender, app_data, user_data) -> None:
        """Открывает форму копирования с безопасными уникальными defaults."""
        mid = str(user_data)
        try:
            name = proj_svc.suggest_duplicate_name(mid, self._state.db_path)
            model_id = proj_svc.suggest_model_id(name, self._state.db_path)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Не удалось подготовить копирование модели %s", mid)
            self._set_status(f"Duplicate failed: {exc}")
            return

        if self._duplicate_win and dpg.does_item_exist(self._duplicate_win):
            dpg.delete_item(self._duplicate_win)
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        with dpg.window(label="Duplicate model", modal=True, no_resize=True,
                        width=540, height=260,
                        pos=(max(0, w // 2 - 270), max(0, h // 2 - 130))) as win:
            self._duplicate_win = self._track_modal(win)
            dpg.add_text(f"Create an independent copy of '{mid}'.")
            dpg.add_text("Stored calculations and GUI tabs are not copied.")
            source = self._state.models.get(mid)
            source_dirty = source is not None and source.dirty
            if source_dirty:
                warning = dpg.add_text(
                    "Source has unsaved changes. The next action will save "
                    "them before creating the copy.")
                dpg.bind_item_theme(warning, self._theme_stale())
            dpg.add_spacer(height=6)
            self._duplicate_ids = {
                "name": dpg.add_input_text(label="Model name", width=430,
                                            default_value=name),
                "id": dpg.add_input_text(label="Model id", width=430,
                                          default_value=model_id),
            }
            dpg.add_spacer(height=6)
            with dpg.group(horizontal=True):
                dpg.add_button(
                    label=("Save source & duplicate" if source_dirty
                           else "Duplicate"),
                    width=190, user_data=(mid, win, source_dirty),
                    callback=self._on_duplicate_model_do)
                dpg.add_button(label="Cancel", width=130,
                               callback=lambda: self._close_duplicate_modal(win))

    def _close_duplicate_modal(self, win: int | None = None) -> None:
        target = win or self._duplicate_win
        if target is not None and dpg.does_item_exist(target):
            dpg.delete_item(target)
        if target == self._duplicate_win:
            self._duplicate_win = None
            self._duplicate_ids = {}

    def _on_duplicate_model_do(self, sender, app_data, user_data) -> None:
        mid, win, save_source = user_data
        source = self._state.models.get(mid)
        if save_source and source is not None and source.dirty:
            try:
                self._state.save_model(mid)
            except Exception as exc:  # noqa: BLE001
                logger.exception("Не удалось сохранить модель перед копированием %s", mid)
                self._set_status(f"Save before duplicate failed: {exc}")
                return
        source = self._state.models.get(mid)
        if source is not None and source.dirty:
            self._set_status(
                f"Save model '{mid}' before duplicating to include unsaved edits.")
            return
        ids = self._duplicate_ids
        try:
            name = str(dpg.get_value(ids["name"])).strip()
            model_id = str(dpg.get_value(ids["id"])).strip()
            new_id = proj_svc.duplicate_model(
                self._state.db_path, mid, new_id=model_id, new_name=name)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Не удалось скопировать модель %s", mid)
            self._set_status(f"Duplicate failed: {exc}")
            return
        self._close_duplicate_modal(win)
        self._selected_project = new_id
        self._state.refresh_model_list()
        self._set_status(
            f"Model '{new_id}' created. Double-click it to open the copy.")

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
        if self._jobs.busy:
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

    # --- импорт E300 --------------------------------------------------------

    def _on_import_e300(self, sender=None, app_data=None, user_data=None) -> None:
        """Файловый диалог выбора E300-дека (.inc)."""
        dialog = dpg.add_file_dialog(
            label="Select Eclipse 300 include file",
            directory_selector=False, show=True, modal=True,
            width=720, height=440,
            callback=self._on_e300_file_chosen,
            cancel_callback=lambda s, a: dpg.delete_item(s))
        self._track_modal(dialog)
        dpg.add_file_extension(".inc", parent=dialog)
        dpg.add_file_extension(".*", parent=dialog)

    def _on_e300_file_chosen(self, sender, app_data) -> None:
        if sender and dpg.does_item_exist(sender):
            dpg.delete_item(sender)
        path = (app_data or {}).get("file_path_name") or ""
        if not path:
            return
        try:
            preview = proj_svc.preview_e300(path)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Предпросмотр E300 не удался")
            self._set_status(f"E300 preview failed: {exc}")
            return
        self._e300_path = path
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        with dpg.window(label="Import Eclipse 300", modal=True, no_resize=True,
                        width=620, height=390,
                        pos=(max(0, w // 2 - 310), max(0, h // 2 - 195))) as win:
            self._e300_win = self._track_modal(win)
            dpg.add_text(f"File: {Path(path).name}")
            dpg.add_text(f"Deck EOS: {preview['deck_eos'] or '-'}")
            dpg.add_text(f"Recognized components: {len(preview['recognized'])}")
            dpg.add_text(f"Composition sum: {preview['total_zi']:.8g}")
            if preview["unrecognized"]:
                warn = dpg.add_text(
                    "Unrecognized components will be skipped:\n"
                    + ", ".join(preview["unrecognized"])
                    + f"\nSkipped mole fraction: {preview['skipped_zi']:.8g}",
                    wrap=580)
                dpg.bind_item_theme(warn, self._theme_stale())
            for warning in preview["warnings"]:
                item = dpg.add_text("Warning: " + warning, wrap=580)
                dpg.bind_item_theme(item, self._theme_stale())
            dpg.add_spacer(height=8)
            dpg.add_text("Import creates a new model. Review warnings before continuing.")
            with dpg.group(horizontal=True):
                dpg.add_button(label="Import", width=140,
                               callback=self._on_e300_import_confirm)
                dpg.add_button(label="Cancel", width=120,
                               callback=lambda: dpg.delete_item(win))

    def _on_e300_import_confirm(self, sender=None, app_data=None, user_data=None) -> None:
        path = self._e300_path
        if not path:
            return
        try:
            res = proj_svc.import_e300(
                self._state.db_path, path, allow_unrecognized=True)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Импорт E300 не удался")
            self._set_status(f"E300 import failed: {exc}")
            return
        if self._e300_win and dpg.does_item_exist(self._e300_win):
            dpg.delete_item(self._e300_win)
        self._e300_win = None
        self._e300_path = None
        self._state.refresh_model_list()
        mid = res["model_id"]
        self._restored_models.add(mid)     # свежая модель, восстанавливать нечего
        self._expanded_models.add(mid)
        self._state.enter_model(mid)
        n_bad = len(res["unrecognized"])
        note = (f", {n_bad} unrecognized skipped ({', '.join(res['unrecognized'][:4])})"
                if n_bad else "")
        self._set_status(
            f"Imported E300 deck '{res['name']}' "
            f"({len(res['recognized'])} components, deck EOS {res['deck_eos']})" + note)

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
        path = self._excel_path
        if grp is None or path is None or not dpg.does_item_exist(grp):
            return
        dpg.delete_item(grp, children_only=True)
        sheet = dpg.get_value(self._excel_sheet_id)
        header = bool(dpg.get_value(self._excel_header_id))
        try:
            prev = proj_svc.excel_preview(path, sheet, header)
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
        path = self._excel_path
        if path is None:
            self._set_status("Excel import path is missing.")
            return
        sheet = dpg.get_value(self._excel_sheet_id)
        header = bool(dpg.get_value(self._excel_header_id))
        try:
            res = proj_svc.import_excel(path, header, sheet)
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
            res["recognized"], name=Path(path).stem,
            unrecognized=res["unrecognized"])
        # Открыть форму на СЛЕДУЮЩЕМ кадре: новая modal-модалка, созданная в
        # том же кадре, где закрыли предыдущую (Excel-предпросмотр), у ImGui
        # остаётся скрытой (OpenPopup не срабатывает, пока закрывается старый
        # popup) — из-за этого форма «не появлялась» после загрузки из листа.
        dpg.set_frame_callback(dpg.get_frame_count() + 1,
                               lambda *_: self._open_new_fluid_modal())
        n_bad = len(res["unrecognized"])
        self._set_status(
            f"Imported {len(res['recognized'])} components from Excel"
            + (f", {n_bad} unrecognized skipped." if n_bad else "."))
