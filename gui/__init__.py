"""
GUI-слой PVTcalc (нативный DearPyGui).

Архитектура: тонкий View (DearPyGui) поверх фреймворк-независимого слоя
состояния (`gui.app_state`) и сервисов (`gui.services`), которые знают
только про `calc_core`. Замена стека (PySide / веб) переписывает только
подпакет `gui.view`, не трогая state/services. См. docs/GUI.md.
"""
