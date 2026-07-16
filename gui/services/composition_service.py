"""
Операции над корневым узлом «Состав» — фреймворк-независимо.

Собирает в одном месте: (1) списки допустимых значений для UI (типы EOS и
методы корреляций свойств C7+ — списки взяты из `calc_core.EOS.BaseEOS`
и `calc_core.PlusComponents.*.get_correlation`); (2) тонкие обёртки над
методами `Composition` для правки состава/свойств/BIP, смены EOS и полного
пересчёта `composition_data`.

Не импортирует DearPyGui — используется из `gui.app_state`, а тесты гоняют
его напрямую (см. tests/test_gui_composition_service.py).
"""

import logging

from calc_core.EOS.BaseEOS import EOSType
from calc_core.Composition.Composition import Composition

logger = logging.getLogger(__name__)

# Допустимые типы EOS (значения `EOSType`).
EOS_OPTIONS: list[str] = [e.value for e in EOSType]

# Методы корреляций C7+ по свойству — источник истины: `get_correlation`
# соответствующего класса в calc_core/PlusComponents/.
CORRELATION_OPTIONS: dict[str, list[str]] = {
    "critical_temperature": [
        "roess", "nokey", "cavett", "kesler lee", "pedersen", "standing",
        "sim daubert", "riazi daubert", "mogoulas tassios", "twu",
        "watansiri owens starling",
    ],
    "critical_pressure": [
        "kesler lee", "riazi daubert", "cavett", "pedersen", "standing",
        "sim daubert", "mogoulas tassios",
    ],
    "acentric_factor": [
        "kesler lee", "riazi al sahhaf", "edmister", "mogoulas tassios",
    ],
    "critical_volume": [
        "hall yarborough", "riazi daubert", "reid", "lohrenz",
    ],
    "Kw": [
        "k watson", "riazi daubert",
    ],
    "shift_parameter": [
        "jhaveri youngren", "srk", "pr",
    ],
}

# Порядок и подписи свойств-корреляций в UI.
CORRELATION_LABELS: list[tuple[str, str]] = [
    ("Kw", "Kw (Watson)"),
    ("critical_temperature", "Tc"),
    ("critical_pressure", "Pc"),
    ("acentric_factor", "acentric omega"),
    ("critical_volume", "Vc"),
    ("shift_parameter", "shift"),
]

# Дефолтная конфигурация корреляций (совпадает с примером в докстринге
# `Composition.evaluate_composition_data`; имена методов — как в
# `get_correlation`, со строчными буквами и пробелами).
DEFAULT_C7_CORRELATIONS: dict[str, str] = {
    "critical_temperature": "pedersen",
    "critical_pressure": "riazi daubert",
    "acentric_factor": "riazi al sahhaf",
    "critical_volume": "hall yarborough",
    "Kw": "k watson",
    "shift_parameter": "jhaveri youngren",
}

# Свойства компонент, доступные для точечной правки (подмножество
# `Composition.edit_component_properties`), в том же порядке, что колонки
# таблицы состава во View.
EDITABLE_PROPERTIES: list[str] = [
    "molar_mass", "critical_temperature", "critical_pressure",
    "acentric_factor", "critical_volume", "shift_parameter", "Kw",
]


def has_c7_plus(composition: Composition) -> bool:
    """True, если в составе есть хотя бы одна C7+ компонента (нужны корреляции)."""
    flags = composition.composition_data.get("c7_plus_flag", {})
    return any(bool(v) for v in flags.values())


def default_correlations() -> dict[str, str]:
    """Свежая копия дефолтной конфигурации корреляций."""
    return dict(DEFAULT_C7_CORRELATIONS)


def set_eos(composition: Composition, eos_value: str) -> None:
    """
    Меняет тип EOS состава (`Composition.eos_name` setter сам пересчитывает
    EOS-зависимые свойства и BIP). Значение — одно из `EOS_OPTIONS`.
    """
    composition.eos_name = EOSType(eos_value)
    logger.info("EOS изменён на %s", eos_value)


def recalculate(composition: Composition, correlations: dict[str, str]) -> None:
    """
    Полный пересчёт `composition_data` под выбранные корреляции C7+
    (`Composition.evaluate_composition_data`). EOS должен быть уже выставлен.
    """
    composition.evaluate_composition_data(c7_plus_correlations=correlations)
    logger.info("Свойства состава пересчитаны (корреляции: %s)", correlations)


def edit_zi(composition: Composition, component: str, value: float) -> None:
    """Меняет мольную долю компонента (без нормализации — см. `normalize`)."""
    composition.composition[component] = float(value)


def normalize(composition: Composition) -> None:
    """Нормирует мольные доли так, чтобы их сумма была равна 1."""
    composition.normalize_composition()


def sum_zi(composition: Composition) -> float:
    """Текущая сумма мольных долей (для индикации ненормированности)."""
    return sum(composition.composition.values())


def edit_property(composition: Composition, component: str,
                  property_name: str, value: float) -> None:
    """
    Точечно правит свойство компонента (`edit_component_properties`); для
    Tc/Pc/omega движок сам пересчитывает EOS-зависимые параметры этой
    компоненты.
    """
    composition.edit_component_properties(component, property_name, float(value))


def edit_bip(composition: Composition, comp1: str, comp2: str, value: float) -> None:
    """Симметрично задаёт коэффициент парного взаимодействия пары компонент."""
    composition.edit_bip_for_components(comp1, comp2, float(value))


def get_bip(composition: Composition, comp1: str, comp2: str) -> float:
    """Текущее значение BIP для пары (0.0, если не задано)."""
    return composition.composition_data.get("bip", {}).get(comp1, {}).get(comp2, 0.0)
