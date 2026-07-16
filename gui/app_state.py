"""
Фреймворк-независимое состояние приложения (ViewModel).

Этот модуль **не импортирует DearPyGui** — только `calc_core` и сервисы
`gui.services`. Здесь описан вычислительный граф `Проект → Модель →
Варианты → Узлы` и статусы узлов (актуален / устарел / считается), на
который опирается реактивная стратегия «ленивая инвалидация + пересчитать
всё» (см. docs/GUI.md). View читает это состояние и вызывает его команды.

В Фазе 0+1 наполнены загрузка списка моделей из `models.json` и открытие
модели (корневой узел «Состав»). Узлы экспериментов (флэш и т.д.) заведены
как заготовки — их вычисление добавляется в последующих фазах.
"""

import logging
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional

from calc_core.Composition.Composition import Composition

from gui.services import composition_service as comp_svc
from gui.services.model_repository import ModelRepository, ModelSummary

logger = logging.getLogger(__name__)


class NodeStatus(Enum):
    """Статус узла вычислительного графа в терминах реактивности."""

    EMPTY = auto()      # ещё не считался
    FRESH = auto()      # результат актуален
    STALE = auto()      # вход изменился, требуется пересчёт
    RUNNING = auto()    # идёт вычисление


class NodeKind(Enum):
    """Тип узла графа."""

    COMPOSITION = auto()   # корневой узел: состав + свойства
    FLASH = auto()         # флэш в точке P,T
    # CCE / DLE / PHASE_ENVELOPE / ... — добавляются в следующих фазах


@dataclass
class GraphNode:
    """
    Узел вычислительного графа модели.

    `result` намеренно нетипизирован (`object`) — конкретный тип зависит от
    вида узла (`FlashResult` для флэша и т.д.); state-слой не привязан к
    структурам результатов движка.
    """

    node_id: str
    kind: NodeKind
    title: str
    status: NodeStatus = NodeStatus.EMPTY
    params: dict = field(default_factory=dict)
    result: object = None
    error: Optional[str] = None
    # id узлов-предков (правка предка помечает потомков STALE)
    upstream: list[str] = field(default_factory=list)


@dataclass
class Variant:
    """
    Вариант модели = форк с другими свойствами / корреляциями / EOS.

    В Фазе 0+1 у каждой модели заводится единственный вариант `Base` с
    корневым узлом состава; наследование вариантов — открытый вопрос
    следующих сессий (см. docs/GUI.md).
    """

    variant_id: str
    title: str
    composition: Optional[Composition] = None
    nodes: dict[str, GraphNode] = field(default_factory=dict)


@dataclass
class Model:
    """Модель флюида: метаданные из `models.json` + загруженные варианты."""

    model_id: str
    title: str
    field_name: Optional[str] = None
    eos: Optional[str] = None
    n_components: int = 0
    loaded: bool = False
    variants: dict[str, Variant] = field(default_factory=dict)


class AppState:
    """
    Корень состояния приложения. Держит список моделей проекта и активный
    выбор (модель / вариант / узел). View подписан на изменения через
    коллбэки `on_change` (простой наблюдатель без внешних зависимостей).
    """

    def __init__(self, repository: ModelRepository):
        self._repo = repository
        self.models: dict[str, Model] = {}
        self.active_model_id: Optional[str] = None
        self.active_variant_id: Optional[str] = None
        self.active_node_id: Optional[str] = None
        self._listeners: list = []

    # --- наблюдатель -----------------------------------------------------

    def subscribe(self, callback) -> None:
        """Регистрирует коллбэк без аргументов, вызываемый при изменении состояния."""
        self._listeners.append(callback)

    def _notify(self) -> None:
        for cb in self._listeners:
            cb()

    # --- загрузка списка моделей ----------------------------------------

    def refresh_model_list(self) -> None:
        """Читает сводку моделей из репозитория (без загрузки составов)."""
        summaries: list[ModelSummary] = self._repo.list_models()
        # сохраняем уже загруженные варианты, если модель осталась в файле
        existing = self.models
        self.models = {}
        for s in summaries:
            prev = existing.get(s.model_id)
            model = Model(
                model_id=s.model_id,
                title=s.title,
                field_name=s.field_name,
                eos=s.eos,
                n_components=s.n_components,
            )
            if prev is not None and prev.loaded:
                model.loaded = True
                model.variants = prev.variants
            self.models[s.model_id] = model
        logger.info("Список моделей обновлён: %d шт.", len(self.models))
        self._notify()

    # --- открытие модели -------------------------------------------------

    def open_model(self, model_id: str) -> None:
        """
        Загружает состав модели из `models.json` (лениво, один раз) и делает
        её активной, создавая вариант `Base` с корневым узлом «Состав».
        """
        model = self.models.get(model_id)
        if model is None:
            raise KeyError(f"Модель '{model_id}' отсутствует в списке")

        if not model.loaded:
            composition = self._repo.load_composition(model_id)
            base = Variant(variant_id="base", title="Base", composition=composition)
            base.nodes["composition"] = GraphNode(
                node_id="composition",
                kind=NodeKind.COMPOSITION,
                title="Composition",
                status=NodeStatus.FRESH,
                params={
                    "eos": composition.eos_name.value,
                    "correlations": comp_svc.default_correlations(),
                },
            )
            model.variants["base"] = base
            model.loaded = True
            logger.info("Модель '%s' загружена (%d компонентов)",
                        model_id, model.n_components)

        self.active_model_id = model_id
        self.active_variant_id = "base"
        self.active_node_id = "composition"
        self._notify()

    # --- команды корневого узла «Состав» --------------------------------
    #
    # Immediate-правки (zi/свойство/BIP применяются движком сразу и остаются
    # консистентными) НЕ дёргают наблюдателя — иначе View пересобрал бы поля
    # ввода и потерял фокус. Дискретные действия (EOS/корреляции/пересчёт/
    # нормировка) вызывают `_notify()` и полную перерисовку.

    @property
    def composition_node(self) -> Optional[GraphNode]:
        """Корневой узел «Состав» активного варианта (или None)."""
        variant = self.active_variant
        if variant is None:
            return None
        return variant.nodes.get("composition")

    def set_composition_eos(self, eos_value: str) -> None:
        """Меняет EOS активного состава (применяется сразу) и перерисовывает."""
        composition = self.active_composition
        node = self.composition_node
        if composition is None or node is None:
            return
        comp_svc.set_eos(composition, eos_value)
        node.params["eos"] = eos_value
        self._notify()

    def set_correlation(self, property_name: str, method: str) -> None:
        """
        Запоминает выбор корреляции C7+ и помечает узел `STALE` — значение
        вступит в силу после `recalculate_composition()`.
        """
        node = self.composition_node
        if node is None:
            return
        node.params.setdefault("correlations", comp_svc.default_correlations())
        node.params["correlations"][property_name] = method
        node.status = NodeStatus.STALE
        self._notify()

    def recalculate_composition(self) -> None:
        """Полный пересчёт свойств под выбранные EOS+корреляции; узел → `FRESH`."""
        composition = self.active_composition
        node = self.composition_node
        if composition is None or node is None:
            return
        node.status = NodeStatus.RUNNING
        self._notify()
        try:
            comp_svc.set_eos(composition, node.params.get("eos", composition.eos_name.value))
            comp_svc.recalculate(composition, node.params.get(
                "correlations", comp_svc.default_correlations()))
            node.status = NodeStatus.FRESH
            node.error = None
        except Exception as exc:  # noqa: BLE001 — донести ошибку до UI
            node.status = NodeStatus.STALE
            node.error = str(exc)
            logger.exception("Пересчёт свойств состава не удался")
        self._notify()

    def edit_zi(self, component: str, value: float) -> None:
        """Immediate-правка мольной доли (без перерисовки)."""
        composition = self.active_composition
        if composition is not None:
            comp_svc.edit_zi(composition, component, value)

    def normalize_composition(self) -> None:
        """Нормировка мольных долей + перерисовка (значения в полях меняются)."""
        composition = self.active_composition
        if composition is not None:
            comp_svc.normalize(composition)
            self._notify()

    def edit_component_property(self, component: str, property_name: str,
                               value: float) -> None:
        """Immediate-правка свойства компонента (без перерисовки)."""
        composition = self.active_composition
        if composition is not None:
            comp_svc.edit_property(composition, component, property_name, value)

    def edit_bip(self, comp1: str, comp2: str, value: float) -> None:
        """Immediate-правка BIP пары компонент (без перерисовки)."""
        composition = self.active_composition
        if composition is not None:
            comp_svc.edit_bip(composition, comp1, comp2, value)

    # --- удобные геттеры активного выбора -------------------------------

    @property
    def active_model(self) -> Optional[Model]:
        if self.active_model_id is None:
            return None
        return self.models.get(self.active_model_id)

    @property
    def active_variant(self) -> Optional[Variant]:
        model = self.active_model
        if model is None or self.active_variant_id is None:
            return None
        return model.variants.get(self.active_variant_id)

    @property
    def active_composition(self) -> Optional[Composition]:
        variant = self.active_variant
        return None if variant is None else variant.composition
