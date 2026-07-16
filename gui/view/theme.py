"""
Светлая глобальная тема приложения.

В DearPyGui нет готового «light»-переключателя (по умолчанию тёмная тема),
поэтому палитра задаётся явно — по мотивам Dear ImGui `StyleColorsLight`.
Цвета сопоставляются с константами `mvThemeCol_*` через `getattr` с
пропуском отсутствующих (имена некоторых констант различаются между
версиями DPG — так тема не падает на незнакомом имени).

Использование: после `create_context()` вызвать
`dpg.bind_theme(build_light_theme())`.
"""

import dearpygui.dearpygui as dpg

# {имя константы mvThemeCol_*: (R, G, B, A)} в диапазоне 0..255.
# Значения — Dear ImGui StyleColorsLight, переведённые из 0..1.
_LIGHT_COLORS: dict[str, tuple[int, int, int, int]] = {
    "mvThemeCol_Text": (0, 0, 0, 255),
    "mvThemeCol_TextDisabled": (153, 153, 153, 255),
    "mvThemeCol_WindowBg": (240, 240, 240, 255),
    "mvThemeCol_ChildBg": (0, 0, 0, 0),
    "mvThemeCol_PopupBg": (255, 255, 255, 250),
    "mvThemeCol_Border": (0, 0, 0, 77),
    "mvThemeCol_BorderShadow": (0, 0, 0, 0),
    "mvThemeCol_FrameBg": (255, 255, 255, 255),
    "mvThemeCol_FrameBgHovered": (66, 150, 250, 102),
    "mvThemeCol_FrameBgActive": (66, 150, 250, 171),
    "mvThemeCol_TitleBg": (245, 245, 245, 255),
    "mvThemeCol_TitleBgActive": (209, 209, 209, 255),
    "mvThemeCol_TitleBgCollapsed": (255, 255, 255, 130),
    "mvThemeCol_MenuBarBg": (219, 219, 219, 255),
    "mvThemeCol_ScrollbarBg": (250, 250, 250, 135),
    "mvThemeCol_ScrollbarGrab": (176, 176, 176, 204),
    "mvThemeCol_ScrollbarGrabHovered": (125, 125, 125, 204),
    "mvThemeCol_ScrollbarGrabActive": (125, 125, 125, 255),
    "mvThemeCol_CheckMark": (66, 150, 250, 255),
    "mvThemeCol_SliderGrab": (66, 150, 250, 199),
    "mvThemeCol_SliderGrabActive": (117, 138, 204, 153),
    "mvThemeCol_Button": (66, 150, 250, 102),
    "mvThemeCol_ButtonHovered": (66, 150, 250, 255),
    "mvThemeCol_ButtonActive": (15, 135, 250, 255),
    "mvThemeCol_Header": (66, 150, 250, 79),
    "mvThemeCol_HeaderHovered": (66, 150, 250, 204),
    "mvThemeCol_HeaderActive": (66, 150, 250, 255),
    "mvThemeCol_Separator": (99, 99, 99, 158),
    "mvThemeCol_SeparatorHovered": (36, 112, 204, 199),
    "mvThemeCol_SeparatorActive": (36, 112, 204, 255),
    "mvThemeCol_ResizeGrip": (89, 89, 89, 43),
    "mvThemeCol_ResizeGripHovered": (66, 150, 250, 171),
    "mvThemeCol_ResizeGripActive": (66, 150, 250, 242),
    "mvThemeCol_Tab": (194, 204, 214, 237),
    "mvThemeCol_TabHovered": (66, 150, 250, 204),
    "mvThemeCol_TabActive": (153, 186, 224, 255),
    "mvThemeCol_TabUnfocused": (235, 237, 240, 252),
    "mvThemeCol_TabUnfocusedActive": (189, 209, 232, 255),
    "mvThemeCol_PlotLines": (99, 99, 99, 255),
    "mvThemeCol_PlotLinesHovered": (255, 110, 89, 255),
    "mvThemeCol_PlotHistogram": (230, 179, 0, 255),
    "mvThemeCol_PlotHistogramHovered": (255, 115, 0, 255),
    "mvThemeCol_TableHeaderBg": (199, 222, 250, 255),
    "mvThemeCol_TableBorderStrong": (145, 145, 163, 255),
    "mvThemeCol_TableBorderLight": (173, 173, 189, 255),
    "mvThemeCol_TableRowBg": (0, 0, 0, 0),
    "mvThemeCol_TableRowBgAlt": (77, 77, 77, 23),
    "mvThemeCol_TextSelectedBg": (66, 150, 250, 89),
    "mvThemeCol_DragDropTarget": (66, 150, 250, 242),
    "mvThemeCol_NavHighlight": (66, 150, 250, 204),
    "mvThemeCol_NavWindowingHighlight": (179, 179, 179, 179),
    "mvThemeCol_NavWindowingDimBg": (51, 51, 51, 51),
    "mvThemeCol_ModalWindowDimBg": (51, 51, 51, 89),
}


def build_light_theme() -> int:
    """
    Создаёт и возвращает id глобальной светлой темы (применять через
    `dpg.bind_theme(...)`). Должна вызываться после `create_context()`.
    """
    with dpg.theme() as theme_id:
        with dpg.theme_component(dpg.mvAll):
            for const_name, rgba in _LIGHT_COLORS.items():
                const = getattr(dpg, const_name, None)
                if const is not None:
                    dpg.add_theme_color(const, rgba, category=dpg.mvThemeCat_Core)
    return theme_id
