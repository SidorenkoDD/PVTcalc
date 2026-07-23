"""
Светлая глобальная тема приложения.

В DearPyGui нет готового «light»-переключателя (по умолчанию тёмная тема),
поэтому палитра задаётся явно. За основу взята спокойная серо-голубая
гамма таблицы LAB DATA: она сохраняет светлую рабочую поверхность и не
создаёт резких акцентов в формах, дереве и таблицах.
Цвета сопоставляются с константами `mvThemeCol_*` через `getattr` с
пропуском отсутствующих (имена некоторых констант различаются между
версиями DPG — так тема не падает на незнакомом имени).

Использование: после `create_context()` вызвать
`dpg.bind_theme(build_light_theme())`.
"""

import dearpygui.dearpygui as dpg

# {имя константы mvThemeCol_*: (R, G, B, A)} в диапазоне 0..255.
# Значения образуют общую спокойную серо-голубую палитру приложения.
_LIGHT_COLORS: dict[str, tuple[int, int, int, int]] = {
    "mvThemeCol_Text": (35, 45, 55, 255),
    "mvThemeCol_TextDisabled": (120, 133, 145, 255),
    # Нейтральная рабочая поверхность: цвет фона оставлен прежним, чтобы
    # серо-голубые поля и выделения читались как отдельные элементы.
    "mvThemeCol_WindowBg": (240, 240, 240, 255),
    "mvThemeCol_ChildBg": (0, 0, 0, 0),
    "mvThemeCol_PopupBg": (238, 243, 248, 255),
    "mvThemeCol_Border": (143, 160, 177, 180),
    "mvThemeCol_BorderShadow": (0, 0, 0, 0),
    "mvThemeCol_FrameBg": (222, 229, 236, 255),
    "mvThemeCol_FrameBgHovered": (213, 225, 237, 255),
    "mvThemeCol_FrameBgActive": (205, 219, 232, 255),
    "mvThemeCol_TitleBg": (218, 226, 234, 255),
    "mvThemeCol_TitleBgActive": (194, 207, 219, 255),
    "mvThemeCol_TitleBgCollapsed": (222, 229, 236, 180),
    "mvThemeCol_MenuBarBg": (213, 221, 230, 255),
    "mvThemeCol_ScrollbarBg": (226, 233, 240, 255),
    "mvThemeCol_ScrollbarGrab": (177, 196, 214, 255),
    "mvThemeCol_ScrollbarGrabHovered": (143, 160, 177, 255),
    "mvThemeCol_ScrollbarGrabActive": (112, 145, 177, 255),
    "mvThemeCol_CheckMark": (112, 145, 177, 255),
    "mvThemeCol_SliderGrab": (143, 160, 177, 255),
    "mvThemeCol_SliderGrabActive": (112, 145, 177, 255),
    "mvThemeCol_Button": (222, 229, 236, 255),
    "mvThemeCol_ButtonHovered": (213, 225, 237, 255),
    "mvThemeCol_ButtonActive": (195, 212, 228, 255),
    "mvThemeCol_Header": (222, 229, 236, 190),
    "mvThemeCol_HeaderHovered": (213, 225, 237, 255),
    "mvThemeCol_HeaderActive": (195, 212, 228, 255),
    "mvThemeCol_Separator": (177, 196, 214, 255),
    "mvThemeCol_SeparatorHovered": (143, 160, 177, 255),
    "mvThemeCol_SeparatorActive": (112, 145, 177, 255),
    "mvThemeCol_ResizeGrip": (177, 196, 214, 120),
    "mvThemeCol_ResizeGripHovered": (143, 160, 177, 190),
    "mvThemeCol_ResizeGripActive": (112, 145, 177, 255),
    "mvThemeCol_Tab": (214, 223, 231, 255),
    "mvThemeCol_TabHovered": (205, 219, 232, 255),
    "mvThemeCol_TabActive": (194, 207, 219, 255),
    "mvThemeCol_TabUnfocused": (226, 233, 240, 255),
    "mvThemeCol_TabUnfocusedActive": (205, 219, 232, 255),
    "mvThemeCol_PlotLines": (82, 98, 115, 255),
    "mvThemeCol_PlotLinesHovered": (255, 110, 89, 255),
    "mvThemeCol_PlotHistogram": (230, 179, 0, 255),
    "mvThemeCol_PlotHistogramHovered": (255, 115, 0, 255),
    "mvThemeCol_TableHeaderBg": (194, 207, 219, 255),
    "mvThemeCol_TableBorderStrong": (143, 160, 177, 255),
    "mvThemeCol_TableBorderLight": (190, 202, 214, 255),
    "mvThemeCol_TableRowBg": (228, 234, 240, 255),
    "mvThemeCol_TableRowBgAlt": (219, 227, 235, 255),
    "mvThemeCol_TextSelectedBg": (195, 212, 228, 190),
    "mvThemeCol_DragDropTarget": (112, 145, 177, 255),
    "mvThemeCol_NavHighlight": (143, 160, 177, 220),
    "mvThemeCol_NavWindowingHighlight": (177, 196, 214, 200),
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
        # Поля ввода остаются отдельными рабочими поверхностями, но используют
        # ту же палитру, что кнопки и ячейки таблицы LAB DATA.
        with dpg.theme_component(dpg.mvInputText):
            dpg.add_theme_color(dpg.mvThemeCol_FrameBg, (222, 229, 236, 255))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgHovered,
                                (213, 225, 237, 255))
            dpg.add_theme_color(dpg.mvThemeCol_FrameBgActive,
                                (205, 219, 232, 255))
            dpg.add_theme_color(dpg.mvThemeCol_Border, (143, 160, 177, 255))
            dpg.add_theme_style(dpg.mvStyleVar_FrameBorderSize, 1)
    return theme_id
