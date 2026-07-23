"""Общие модальные диалоги экспорта и настроек GUI."""

import logging

import dearpygui.dearpygui as dpg

from gui.services import export_service as exp_out_svc
from gui.services import model_report_service as report_svc
from gui.services import settings_service as settings_svc
from gui.view.contracts import ContextBoundView

logger = logging.getLogger(__name__)


class DialogsViewMixin(ContextBoundView):
    """Диалоги, доступные из главного меню приложения."""

    _selected_project: str | None
    _export_comp: object | None
    _export_label: str
    _export_fmt: str
    _export_eos: str
    _export_ids: dict[str, int]
    _export_win: int | None
    _report_model_id: str | None
    _report_options: dict[str, bool]
    _report_ids: dict[str, int]
    _report_win: int | None
    _settings_ids: dict[str, int]
    _settings_win: int | None

    def _composition_for_export(self):
        """Состав для экспорта: активная модель, иначе выбранная на Projects."""
        if self._state.active_composition is not None:
            model = self._state.active_model
            return self._state.active_composition, (model.title if model else "model")
        mid = self._selected_project
        if mid:
            try:
                return self._state.peek_composition(mid), mid
            except Exception:  # noqa: BLE001
                logger.exception("Не удалось прочитать состав %s для экспорта", mid)
        return None, None

    def _on_open_export_dialog(self, sender=None, app_data=None, user_data=None) -> None:
        comp, label = self._composition_for_export()
        if comp is None:
            self._set_status("Open a model (or select one on Projects) to export.")
            return
        self._export_comp = comp
        self._export_label = label or "model"
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        with dpg.window(label=f"Export - {self._export_label}", modal=True,
                        no_resize=True, no_collapse=True, width=400, height=230,
                        pos=(max(0, w // 2 - 200), max(0, h // 2 - 115))) as win:
            self._track_modal(win)
            self._export_win = win
            ids: dict = {}
            ids["format"] = dpg.add_combo(
                items=exp_out_svc.format_labels(), width=240, label="Format",
                default_value=exp_out_svc.format_labels()[0],
                callback=self._on_export_format_changed)
            dpg.add_separator()
            with dpg.group() as eos_grp:
                ids["eos"] = dpg.add_combo(
                    items=exp_out_svc.FORMATS["e300"]["eos_choices"], width=120,
                    default_value="MPR", label="EOS keyword")
                dpg.add_text("Brusilovskiy has no direct E300 EOS -\n"
                             "pick PR (MPR) or SRK for the deck.")
            ids["_eos_grp"] = eos_grp
            self._export_ids = ids
            dpg.add_spacer(height=6)
            with dpg.group(horizontal=True):
                dpg.add_button(label="Choose file & export", width=180,
                               callback=self._on_export_choose_file)
                dpg.add_button(label="Cancel", width=100,
                               callback=lambda: dpg.delete_item(win))

    def _on_export_format_changed(self, sender, app_data) -> None:
        """Показывает опции EOS только для форматов, где они есть (E300)."""
        fmt = exp_out_svc.format_key_by_label(app_data)
        has_eos = "eos_choices" in exp_out_svc.FORMATS.get(fmt, {})
        if dpg.does_item_exist(self._export_ids.get("_eos_grp", 0)):
            dpg.configure_item(self._export_ids["_eos_grp"], show=has_eos)

    def _on_export_choose_file(self, sender=None, app_data=None, user_data=None) -> None:
        self._export_fmt = exp_out_svc.format_key_by_label(
            dpg.get_value(self._export_ids["format"]))
        self._export_eos = dpg.get_value(self._export_ids.get("eos", 0)) or "MPR"
        ext = exp_out_svc.FORMATS[self._export_fmt]["ext"]
        safe = "".join(c if c.isalnum() or c in "-_" else "_"
                       for c in self._export_label) or "model"
        dialog = dpg.add_file_dialog(
            label=f"Export {self._export_label} to {ext}",
            directory_selector=False, show=True, modal=True,
            width=720, height=440, default_filename=safe,
            callback=self._on_export_file_chosen,
            cancel_callback=lambda s, a: dpg.delete_item(s))
        self._track_modal(dialog)
        dpg.add_file_extension(ext, parent=dialog)
        dpg.add_file_extension(".*", parent=dialog)

    def _on_export_file_chosen(self, sender, app_data) -> None:
        if sender and dpg.does_item_exist(sender):
            dpg.delete_item(sender)
        path = (app_data or {}).get("file_path_name") or ""
        if not path:
            return
        try:
            written = exp_out_svc.export_model(
                self._export_comp, self._export_fmt, path,
                eos_keyword=self._export_eos)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Экспорт не удался")
            self._set_status(f"Export failed: {exc}")
            return
        if self._export_win and dpg.does_item_exist(self._export_win):
            dpg.delete_item(self._export_win)
        self._set_status(f"Exported to {written}")

    # --- инженерный Excel-отчёт ------------------------------------------

    def _on_open_model_report_dialog(self, sender=None, app_data=None, user_data=None) -> None:
        """Opens report selection for the model root selected in the tree."""
        model_id = str(user_data)
        try:
            self._restore_workspace(model_id)
            model = self._state.ensure_model_loaded(model_id)
        except KeyError:
            self._set_status("Model is no longer available.")
            return
        variant = model.variants.get("base")
        if variant is None:
            self._set_status("Model workspace is unavailable.")
            return
        self._report_model_id = model_id
        self._report_options = {
            "composition": True, "flash": True, "dle": True, "cce": True,
            "separator": True, "envelope": True, "lab_data": True,
        }
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        with dpg.window(label=f"Excel report - {model.title}", modal=True,
                        no_resize=True, no_collapse=True, width=430, height=380,
                        pos=(max(0, w // 2 - 215), max(0, h // 2 - 190))) as win:
            self._track_modal(win)
            self._report_win = win
            dpg.add_text("Choose report sections:")
            ids = {}
            ids["composition"] = dpg.add_checkbox(
                label="Composition and component properties", default_value=True)
            dpg.add_separator()
            dpg.add_text("Completed calculations")
            for key, label in (
                ("flash", "Flash"), ("dle", "DLE"), ("cce", "CCE"),
                ("separator", "Separator"), ("envelope", "Phase envelope"),
            ):
                ids[key] = dpg.add_checkbox(label=label, default_value=True)
            dpg.add_separator()
            ids["lab_data"] = dpg.add_checkbox(
                label="Linked laboratory data on experiment sheets", default_value=True)
            self._report_ids = ids
            dpg.add_spacer(height=6)
            with dpg.group(horizontal=True):
                dpg.add_button(label="Choose file & export", width=170,
                               callback=self._on_model_report_prepare_export)
                dpg.add_button(label="Cancel", width=100,
                               callback=lambda: self._close_tracked_modal(win))

    def _on_model_report_prepare_export(self, sender=None, app_data=None,
                                        user_data=None) -> None:
        model_id = self._report_model_id
        model = self._state.models.get(model_id) if model_id else None
        if model is None:
            self._set_status("Model is no longer available.")
            return
        self._report_options = {
            key: bool(dpg.get_value(item_id))
            for key, item_id in self._report_ids.items()
        }
        if self._report_win is not None and dpg.does_item_exist(self._report_win):
            self._close_tracked_modal(self._report_win)
        self._report_win = None
        if model.dirty:
            self._open_model_report_dirty_dialog(model)
            return
        self._open_model_report_file_dialog()

    def _open_model_report_dirty_dialog(self, model) -> None:
        with dpg.window(label="Unsaved model changes", modal=True, no_collapse=True,
                        width=450, height=170) as win:
            self._track_modal(win)
            dpg.add_text(
                "The report can match the saved model or the current temporary state.",
                wrap=410,
            )
            with dpg.group(horizontal=True):
                dpg.add_button(label="Save & export", width=125, user_data=win,
                               callback=self._on_model_report_save_and_export)
                dpg.add_button(label="Export current snapshot", width=155,
                               user_data=win,
                               callback=self._on_model_report_export_snapshot)
                dpg.add_button(label="Cancel", width=80,
                               callback=lambda: self._close_tracked_modal(win))

    def _on_model_report_save_and_export(self, sender, app_data, user_data) -> None:
        model_id = self._report_model_id
        try:
            self._state.save_model(model_id)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Could not save model before report export")
            self._set_status(f"Save failed: {exc}")
            return
        self._close_tracked_modal(user_data)
        self._open_model_report_file_dialog()

    def _on_model_report_export_snapshot(self, sender, app_data, user_data) -> None:
        self._close_tracked_modal(user_data)
        self._open_model_report_file_dialog()

    def _open_model_report_file_dialog(self) -> None:
        model = self._state.models.get(self._report_model_id or "")
        if model is None:
            return
        safe = "".join(char if char.isalnum() or char in "-_" else "_"
                       for char in model.title) or "model_report"
        dialog = dpg.add_file_dialog(
            label=f"Export Excel report - {model.title}", directory_selector=False,
            show=True, modal=True, width=720, height=440,
            default_filename=f"{safe}_report",
            callback=self._on_model_report_file_chosen,
            cancel_callback=lambda s, a: self._close_tracked_modal(s),
        )
        self._track_modal(dialog)
        dpg.add_file_extension(".xlsx", parent=dialog)

    def _on_model_report_file_chosen(self, sender, app_data) -> None:
        self._close_tracked_modal(sender)
        path = (app_data or {}).get("file_path_name") or ""
        model = self._state.models.get(self._report_model_id or "")
        variant = model.variants.get("base") if model else None
        if not path or model is None or variant is None:
            return
        try:
            written = report_svc.export_model_report(
                path, model, variant, self._state.db_path, self._report_options)
        except Exception as exc:  # noqa: BLE001
            logger.exception("Excel report export failed")
            self._set_status(f"Report export failed: {exc}")
            return
        if self._report_win is not None and dpg.does_item_exist(self._report_win):
            self._close_tracked_modal(self._report_win)
        self._report_win = None
        self._set_status(f"Excel report exported to {written}")

    # --- настройки (константы/условия/критерии сходимости) ------------------

    def _on_open_settings(self, sender=None, app_data=None, user_data=None) -> None:
        values = settings_svc.engine_defaults()
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        height = min(max(360, h - 120), 620)
        with dpg.window(label="Engine constants (read only)", modal=True,
                        no_collapse=True,
                        width=460, height=height,
                        pos=(max(0, w // 2 - 230), 60)) as win:
            self._track_modal(win)
            self._settings_win = win
            self._settings_ids = {}
            note = dpg.add_text(
                "Actual values used by the calculation engine. Editing will be\n"
                "available only after an explicit EngineConfig API is added.",
                wrap=430)
            dpg.bind_item_theme(note, self._theme_stale())
            dpg.add_separator()
            with dpg.child_window(height=height - 140, autosize_x=True):
                for grp in settings_svc.SCHEMA:
                    with dpg.collapsing_header(label=grp["group"], default_open=True):
                        for f in grp["fields"]:
                            key = f["key"]
                            self._settings_ids[key] = dpg.add_input_text(
                                label=f["label"], width=180,
                                readonly=True,
                                default_value=self._fmt_setting(values.get(key), f))
            dpg.add_separator()
            dpg.add_button(label="Close", width=100,
                           callback=lambda: dpg.delete_item(win))

    @staticmethod
    def _fmt_setting(value, field) -> str:
        if value is None:
            return ""
        try:
            return field.get("fmt", "%.6g") % float(value)
        except (TypeError, ValueError):
            return str(value)
