"""Проверки создания глобальной темы DearPyGui."""

import dearpygui.dearpygui as dpg

from gui.view.theme import build_light_theme


def test_build_light_theme_creates_theme() -> None:
    """Глобальная серо-голубая тема создаётся в чистом DPG-контексте."""
    dpg.create_context()
    try:
        theme_id = build_light_theme()
        assert dpg.does_item_exist(theme_id)
    finally:
        dpg.destroy_context()
