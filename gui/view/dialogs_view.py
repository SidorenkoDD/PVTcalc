"""Общие модальные диалоги экспорта и настроек GUI."""

import logging
from typing import Any

import dearpygui.dearpygui as dpg

from gui.services import export_service as exp_out_svc
from gui.services import settings_service as settings_svc

logger = logging.getLogger(__name__)


class DialogsViewMixin:
    """Диалоги, доступные из главного меню приложения."""

    _state: Any
    _selected_project: Any
    _export_comp: Any
    _export_label: Any
    _export_fmt: Any
    _export_eos: Any
    _export_ids: Any
    _export_win: Any
    _settings_ids: Any
    _settings_win: Any
    _set_status: Any
    _theme_stale: Any
    _track_modal: Any

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

    # --- настройки (константы/условия/критерии сходимости) ------------------

    def _on_open_settings(self, sender=None, app_data=None, user_data=None) -> None:
        values = settings_svc.load_settings()
        w, h = dpg.get_viewport_width(), dpg.get_viewport_height()
        height = min(max(360, h - 120), 620)
        with dpg.window(label="Settings", modal=True, no_collapse=True,
                        width=460, height=height,
                        pos=(max(0, w // 2 - 230), 60)) as win:
            self._track_modal(win)
            self._settings_win = win
            self._settings_ids = {}
            note = dpg.add_text(
                "Current engine values. Saved to gui_settings.json, but do NOT\n"
                "affect calculations yet (engine wiring is planned).",
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
                                default_value=self._fmt_setting(values.get(key), f))
            dpg.add_separator()
            with dpg.group(horizontal=True):
                dpg.add_button(label="Save", width=110, callback=self._on_settings_save)
                dpg.add_button(label="Reset to defaults", width=150,
                               callback=self._on_settings_reset)
                dpg.add_button(label="Cancel", width=90,
                               callback=lambda: dpg.delete_item(win))

    @staticmethod
    def _fmt_setting(value, field) -> str:
        if value is None:
            return ""
        try:
            return field.get("fmt", "%.6g") % float(value)
        except (TypeError, ValueError):
            return str(value)

    def _on_settings_save(self, sender=None, app_data=None, user_data=None) -> None:
        values = {}
        try:
            for key, wid in self._settings_ids.items():
                values[key] = float(dpg.get_value(wid))
        except (ValueError, TypeError) as exc:
            self._set_status(f"Invalid settings value: {exc}")
            return
        settings_svc.save_settings(values)
        if self._settings_win and dpg.does_item_exist(self._settings_win):
            dpg.delete_item(self._settings_win)
        self._set_status("Settings saved (not yet applied to calculations).")

    def _on_settings_reset(self, sender=None, app_data=None, user_data=None) -> None:
        """Возвращает поля к дефолтам движка (без сохранения)."""
        defaults = settings_svc.engine_defaults()
        fmt_by_key = {f["key"]: f for grp in settings_svc.SCHEMA for f in grp["fields"]}
        for key, wid in self._settings_ids.items():
            if dpg.does_item_exist(wid):
                dpg.set_value(wid, self._fmt_setting(defaults.get(key),
                                                     fmt_by_key.get(key, {})))
        self._set_status("Settings reset to engine defaults (not saved yet).")
