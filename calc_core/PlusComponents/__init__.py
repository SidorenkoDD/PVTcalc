"""
Не обычный package __init__ — не переэкспортирует ничего из подмодулей
(`__all__` нет), просто содержит несколько неиспользуемых импортов. Похоже
на случайно оставленный код верхнего уровня, а не на намеренный re-export
(см. CLAUDE.md). Не трогать без отдельного запроса — не мешает работе пакета.
"""

import math
from typing import Dict, Callable
from calc_core.Utils.Constants import CONSTANT_R