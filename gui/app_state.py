"""
Фреймворк-независимое состояние приложения (ViewModel).

Этот модуль **не импортирует DearPyGui** — только `calc_core` и сервисы
`gui.services`. Здесь описан вычислительный граф `Проект → Модель →
Вариант → Узлы` и статусы узлов (актуален / устарел / считается), на
который опирается реактивная стратегия «ленивая инвалидация» (см. docs/GUI.md).
View читает это состояние и вызывает его команды.

Узлы одного варианта:
- один корневой узел «Composition» (id `composition`) — состав + свойства;
- ноль и более узлов «Flash» (id `flash_1`, `flash_2`, …) — история
  проведённых флэш-расчётов (каждый запуск = отдельный узел).

Вариант хранит, какие узлы **открыты как вкладки** (`open_node_ids`) и какой
из них активен (`active_node_id`) — так каждая модель помнит свой набор
вкладок при переключении.
"""

import copy
import logging
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional

from calc_core.Composition.Composition import Composition
from calc_core.Utils.ModelStore import ModelStoreError
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

    COMPOSITION = auto()     # корневой узел: состав + свойства
    FLASH = auto()           # флэш в точке P,T
    EXPERIMENT = auto()      # PVT-эксперимент (CCE/DLE/Separator), params["kind"]
    COMPARE = auto()         # сравнение нескольких расчётов (params["members"])
    PHASE_ENVELOPE = auto()  # фазовая огибающая P-T + крит. точка + пластовое Psat
    # SATURATION / ... — добавляются позже


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


@dataclass(frozen=True)
class NodeRef:
    """Стабильный адрес узла, не зависящий от активной вкладки/модели."""

    model_id: str
    variant_id: str
    node_id: str


class StateChangeKind(Enum):
    """Минимальная область View, которую нужно обновить после команды."""

    FULL = auto()
    NAVIGATION = auto()
    MODEL_LIST = auto()
    TREE = auto()
    WORKSPACE = auto()
    NODE = auto()


@dataclass(frozen=True)
class StateChange:
    """Событие состояния для гранулярного рендера DearPyGui."""

    kind: StateChangeKind
    node_ref: Optional[NodeRef] = None


@dataclass
class Variant:
    """
    Вариант модели = форк с другими свойствами / корреляциями / EOS.

    Пока у каждой модели единственный вариант `Base`. Наследование вариантов —
    открытый вопрос следующих сессий (см. docs/GUI.md).
    """

    variant_id: str
    title: str
    composition: Optional[Composition] = None
    nodes: dict[str, GraphNode] = field(default_factory=dict)
    flash_seq: int = 0                                  # счётчик id флэш-узлов
    exp_seq: int = 0                                    # счётчик id узлов экспериментов
    env_seq: int = 0                                    # счётчик id узлов фазовой огибающей
    open_node_ids: list[str] = field(default_factory=list)   # вкладки
    active_node_id: Optional[str] = None                # активная вкладка
    # стеки undo/redo (снимки состояния варианта), у каждой модели свои
    undo_stack: list = field(default_factory=list)
    redo_stack: list = field(default_factory=list)

    def flash_runs(self) -> list[GraphNode]:
        """Все флэш-узлы варианта в порядке добавления (история расчётов)."""
        return [n for n in self.nodes.values() if n.kind is NodeKind.FLASH]

    def experiment_runs(self) -> list[GraphNode]:
        """Все узлы-эксперименты варианта в порядке добавления."""
        return [n for n in self.nodes.values() if n.kind is NodeKind.EXPERIMENT]

    def envelope_runs(self) -> list[GraphNode]:
        """Все узлы фазовой огибающей варианта в порядке добавления."""
        return [n for n in self.nodes.values() if n.kind is NodeKind.PHASE_ENVELOPE]


@dataclass
class Model:
    """Модель флюида: метаданные из `models.json` + загруженные варианты."""

    model_id: str
    title: str
    project_id: str = ""
    field_name: Optional[str] = None
    eos: Optional[str] = None
    n_components: int = 0
    loaded: bool = False
    dirty: bool = False
    variants: dict[str, Variant] = field(default_factory=dict)
    # полная сводка из репозитория (created_at/t_res/results_brief) — для Projects
    summary: Optional[ModelSummary] = None


@dataclass(frozen=True)
class Project:
    """Проект первого уровня, объединяющий несколько моделей-вариантов."""

    project_id: str
    title: str
    model_ids: tuple[str, ...] = ()


class AppState:
    """
    Корень состояния приложения. Держит список моделей проекта и активный
    выбор (модель / вариант / узел). View подписан на изменения через
    коллбэки (простой наблюдатель без внешних зависимостей).
    """

    def __init__(self, repository: ModelRepository):
        self._repo = repository
        self.models: dict[str, Model] = {}
        self.projects: dict[str, Project] = {}
        self.model_list_error: str | None = None
        self.active_model_id: Optional[str] = None
        self.active_project_id: Optional[str] = None
        self.active_variant_id: Optional[str] = None
        # текущий экран: "projects" (стартовый) | "workspace"
        # (форма нового флюида — модальное окно поверх Projects, не экран)
        self.current_screen: str = "projects"
        self._listeners: list[Callable[[], None]] = []
        self._change_listeners: list[Callable[[StateChange], None]] = []
        self._dirty_listeners: list[Callable[[str], None]] = []

    @property
    def db_path(self) -> str:
        """Путь к файлу моделей (models.json) — для сервисов создания/импорта."""
        return str(self._repo.db_path)

    # --- наблюдатель -----------------------------------------------------

    def subscribe(self, callback: Callable[[], None]) -> None:
        """Регистрирует коллбэк без аргументов, вызываемый при изменении состояния."""
        self._listeners.append(callback)

    def subscribe_changes(self, callback: Callable[[StateChange], None]) -> None:
        """Регистрирует слушатель гранулярных событий для View."""
        self._change_listeners.append(callback)

    def subscribe_dirty(self, callback) -> None:
        """Подписывает коллбэк ``callback(model_id)`` на несохранённые правки."""
        self._dirty_listeners.append(callback)

    def _notify(self, change: StateChange | None = None) -> None:
        for cb in self._listeners:
            cb()
        event = change or StateChange(StateChangeKind.FULL)
        for cb in self._change_listeners:
            cb(event)

    def notify(self, change: StateChange | None = None) -> None:
        """Публично сообщает View о пакетном изменении состояния."""
        self._notify(change)

    def _mark_dirty(self) -> None:
        model = self.active_model
        if model is None:
            return
        model.dirty = True
        for cb in self._dirty_listeners:
            cb(model.model_id)

    # --- undo / redo (снимки состояния активного варианта) --------------
    #
    # `_push_undo()` вызывается в НАЧАЛЕ каждой откатываемой операции (правки
    # состава и структурные операции над узлами); запуск флэша и навигация по
    # вкладкам НЕ откатываются. Снимок = глубокая копия узлов + вкладок +
    # состава (deepcopy состава ~1 мс, см. проверку).

    _UNDO_CAP = 40

    @staticmethod
    def _snapshot_variant(variant: "Variant") -> dict:
        # Результаты расчётов считаются неизменяемыми и разделяются между
        # снимками. Копируем только редактируемые метаданные узлов — иначе
        # 40 undo-снимков тяжёлой огибающей занимали бы сотни мегабайт.
        nodes = {
            nid: GraphNode(
                node_id=node.node_id,
                kind=node.kind,
                title=node.title,
                status=node.status,
                params=copy.deepcopy(node.params),
                result=node.result,
                error=node.error,
                upstream=list(node.upstream),
            )
            for nid, node in variant.nodes.items()
        }
        return {
            "nodes": nodes,
            "open": list(variant.open_node_ids),
            "active": variant.active_node_id,
            "flash_seq": variant.flash_seq,
            "exp_seq": variant.exp_seq,
            "env_seq": variant.env_seq,
            "composition": copy.deepcopy(variant.composition),
        }

    @staticmethod
    def _restore_variant(variant: "Variant", snap: dict) -> None:
        variant.nodes = snap["nodes"]
        variant.open_node_ids = list(snap["open"])
        variant.active_node_id = snap["active"]
        variant.flash_seq = snap["flash_seq"]
        variant.exp_seq = snap["exp_seq"]
        variant.env_seq = snap["env_seq"]
        variant.composition = snap["composition"]

    def _push_undo(self) -> None:
        variant = self.active_variant
        if variant is None:
            return
        variant.undo_stack.append(self._snapshot_variant(variant))
        if len(variant.undo_stack) > self._UNDO_CAP:
            variant.undo_stack.pop(0)
        variant.redo_stack.clear()

    def can_undo(self) -> bool:
        v = self.active_variant
        return bool(v and v.undo_stack)

    def can_redo(self) -> bool:
        v = self.active_variant
        return bool(v and v.redo_stack)

    def undo(self) -> None:
        v = self.active_variant
        if v is None or not v.undo_stack:
            return
        v.redo_stack.append(self._snapshot_variant(v))
        self._restore_variant(v, v.undo_stack.pop())
        self._mark_dirty()
        self._notify(StateChange(StateChangeKind.WORKSPACE))

    def redo(self) -> None:
        v = self.active_variant
        if v is None or not v.redo_stack:
            return
        v.undo_stack.append(self._snapshot_variant(v))
        self._restore_variant(v, v.redo_stack.pop())
        self._mark_dirty()
        self._notify(StateChange(StateChangeKind.WORKSPACE))

    def clear_history(self) -> None:
        """Сбрасывает стеки undo/redo активного варианта (напр. после restore)."""
        v = self.active_variant
        if v is not None:
            v.undo_stack.clear()
            v.redo_stack.clear()

    # --- загрузка списка моделей ----------------------------------------

    def refresh_model_list(self) -> None:
        """Читает сводку моделей из репозитория (без загрузки составов)."""
        try:
            summaries: list[ModelSummary] = self._repo.list_models()
        except ModelStoreError as exc:
            # Ошибка пользовательского файла не должна прерывать запуск GUI.
            # Оставляем Projects работоспособным и показываем диагностический
            # баннер; запись намеренно не пытаемся "починить" автоматически.
            self.models = {}
            self.projects = {}
            self.active_model_id = None
            self.active_project_id = None
            self.active_variant_id = None
            self.model_list_error = str(exc)
            logger.error("Не удалось обновить список моделей: %s", exc)
            self._notify(StateChange(StateChangeKind.MODEL_LIST))
            return

        self.model_list_error = None
        existing = self.models
        self.models = {}
        grouped_ids: dict[str, list[str]] = {}
        grouped_titles: dict[str, str] = {}
        for s in summaries:
            prev = existing.get(s.model_id)
            model = Model(
                model_id=s.model_id,
                title=s.title,
                project_id=s.project_id or s.model_id,
                field_name=s.field_name,
                eos=s.eos,
                n_components=s.n_components,
                summary=s,
            )
            if prev is not None and prev.loaded:
                model.loaded = True
                model.variants = prev.variants
                model.dirty = prev.dirty
            self.models[s.model_id] = model
            project_id = model.project_id or model.model_id
            grouped_ids.setdefault(project_id, []).append(model.model_id)
            grouped_titles.setdefault(project_id,
                                       s.project_name or model.title)
        self.projects = {
            project_id: Project(project_id=project_id,
                                title=grouped_titles[project_id],
                                model_ids=tuple(model_ids))
            for project_id, model_ids in grouped_ids.items()
        }
        if self.active_model_id not in self.models:
            self.active_model_id = None
            self.active_project_id = None
            self.active_variant_id = None
        elif self.active_model_id is not None:
            self.active_project_id = self.models[self.active_model_id].project_id
        logger.info("Список моделей обновлён: %d шт.", len(self.models))
        self._notify(StateChange(StateChangeKind.MODEL_LIST))

    # --- загрузка / выбор модели ----------------------------------------

    def _ensure_loaded(self, model_id: str) -> Model:
        """Лениво загружает состав модели и заводит вариант `Base` с узлом состава."""
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
                    "correlations": (self._repo.load_correlations(model_id)
                                     or comp_svc.default_correlations()),
                },
            )
            model.variants["base"] = base
            model.loaded = True
            logger.info("Модель '%s' загружена (%d компонентов)",
                        model_id, model.n_components)
        return model

    def peek_composition(self, model_id: str):
        """
        Загружает состав модели «на посмотреть» (read-only), НЕ помещая его в
        `models` и не меняя активную модель — для предпросмотра на Projects.
        """
        return self._repo.load_composition(model_id)

    def set_active_model(self, model_id: str, *, notify: bool = True) -> None:
        """Делает модель активной (лениво загрузив состав). Вкладки не трогает."""
        self._ensure_loaded(model_id)
        self.active_model_id = model_id
        self.active_project_id = self.models[model_id].project_id
        self.active_variant_id = "base"
        if notify:
            self._notify(StateChange(StateChangeKind.WORKSPACE))

    # обратная совместимость со старым именем
    open_model = set_active_model

    # --- навигация между экранами ----------------------------------------

    def show_projects(self) -> None:
        """Переход на стартовый экран Projects."""
        self.current_screen = "projects"
        self._notify(StateChange(StateChangeKind.NAVIGATION))

    def enter_model(self, model_id: str, *, notify: bool = True) -> None:
        """Открывает модель в рабочем пространстве (Projects → workspace)."""
        self._ensure_loaded(model_id)
        self.active_model_id = model_id
        self.active_project_id = self.models[model_id].project_id
        self.active_variant_id = "base"
        self.current_screen = "workspace"
        if notify:
            self._notify(StateChange(StateChangeKind.NAVIGATION))

    # --- вкладки (открытые узлы) ----------------------------------------

    def open_node(self, node_id: str) -> None:
        """Открывает узел как вкладку (или фокусирует уже открытую), делает активным."""
        variant = self.active_variant
        if variant is None or node_id not in variant.nodes:
            return
        if node_id not in variant.open_node_ids:
            variant.open_node_ids.append(node_id)
        variant.active_node_id = node_id
        self._notify(StateChange(StateChangeKind.WORKSPACE))

    def focus_node(self, node_id: str) -> None:
        """Делает узел активным без перерисовки рабочей области (напр. смена вкладки)."""
        variant = self.active_variant
        if variant is not None and node_id in variant.open_node_ids:
            variant.active_node_id = node_id

    def close_node(self, node_id: str, *, notify: bool = True) -> None:
        """Закрывает вкладку узла; активной становится соседняя (или ни одной)."""
        variant = self.active_variant
        if variant is None or node_id not in variant.open_node_ids:
            return
        idx = variant.open_node_ids.index(node_id)
        variant.open_node_ids.remove(node_id)
        if variant.active_node_id == node_id:
            if variant.open_node_ids:
                variant.active_node_id = variant.open_node_ids[min(
                    idx, len(variant.open_node_ids) - 1)]
            else:
                variant.active_node_id = None
        if notify:
            self._notify(StateChange(StateChangeKind.WORKSPACE))

    # --- команды корневого узла «Состав» --------------------------------
    #
    # Immediate-правки (zi/свойство/BIP) НЕ дёргают наблюдателя — иначе View
    # пересобрал бы поля ввода и потерял фокус. Дискретные действия
    # (EOS/корреляции/пересчёт/нормировка) вызывают `_notify()`.

    @property
    def composition_node(self) -> Optional[GraphNode]:
        """Корневой узел «Состав» активного варианта (или None)."""
        variant = self.active_variant
        return None if variant is None else variant.nodes.get("composition")

    def set_composition_eos(self, eos_value: str) -> None:
        """Меняет EOS активного состава (применяется сразу) и перерисовывает."""
        composition = self.active_composition
        node = self.composition_node
        if composition is None or node is None:
            return
        self._push_undo()
        comp_svc.set_eos(composition, eos_value)
        node.params["eos"] = eos_value
        self._invalidate_flash()
        self._mark_dirty()
        self._notify(StateChange(StateChangeKind.WORKSPACE))

    def set_correlation(self, property_name: str, method: str) -> None:
        """Запоминает выбор корреляции C7+ и помечает узел `STALE` (до пересчёта)."""
        node = self.composition_node
        if node is None:
            return
        self._push_undo()
        node.params.setdefault("correlations", comp_svc.default_correlations())
        node.params["correlations"][property_name] = method
        node.status = NodeStatus.STALE
        self._mark_dirty()
        self._notify(StateChange(StateChangeKind.WORKSPACE))

    def recalculate_composition(self) -> None:
        """Полный пересчёт свойств под выбранные EOS+корреляции; узел → `FRESH`."""
        composition = self.active_composition
        node = self.composition_node
        if composition is None or node is None:
            return
        self._push_undo()
        node.status = NodeStatus.RUNNING
        self._notify(StateChange(StateChangeKind.WORKSPACE))
        try:
            comp_svc.set_eos(composition, node.params.get("eos", composition.eos_name.value))
            comp_svc.recalculate(composition, node.params.get(
                "correlations", comp_svc.default_correlations()))
            node.status = NodeStatus.FRESH
            node.error = None
            self._invalidate_flash()
            self._mark_dirty()
        except Exception as exc:  # noqa: BLE001 — донести ошибку до UI
            node.status = NodeStatus.STALE
            node.error = str(exc)
            logger.exception("Пересчёт свойств состава не удался")
        self._notify(StateChange(StateChangeKind.WORKSPACE))

    def edit_zi(self, component: str, value: float) -> None:
        """Immediate-правка мольной доли (без перерисовки)."""
        composition = self.active_composition
        variant = self.active_variant
        if composition is not None and variant is not None:
            self._push_undo()
            try:
                comp_svc.edit_zi(composition, component, value)
            except Exception:
                variant.undo_stack.pop()
                raise
            self._invalidate_flash()
            self._mark_dirty()

    def normalize_composition(self) -> None:
        """Нормировка мольных долей + перерисовка (значения в полях меняются)."""
        composition = self.active_composition
        if composition is not None:
            self._push_undo()
            comp_svc.normalize(composition)
            self._invalidate_flash()
            self._mark_dirty()
            self._notify(StateChange(StateChangeKind.WORKSPACE))

    def edit_component_property(self, component: str, property_name: str,
                               value: float) -> None:
        """Immediate-правка свойства компонента (без перерисовки)."""
        composition = self.active_composition
        variant = self.active_variant
        if composition is not None and variant is not None:
            self._push_undo()
            try:
                comp_svc.edit_property(composition, component, property_name, value)
            except Exception:
                variant.undo_stack.pop()
                raise
            self._invalidate_flash()
            self._mark_dirty()

    def edit_bip(self, comp1: str, comp2: str, value: float) -> None:
        """Immediate-правка BIP пары компонент (без перерисовки)."""
        composition = self.active_composition
        variant = self.active_variant
        if composition is not None and variant is not None:
            self._push_undo()
            try:
                comp_svc.edit_bip(composition, comp1, comp2, value)
            except Exception:
                variant.undo_stack.pop()
                raise
            self._invalidate_flash()
            self._mark_dirty()

    def save_model(self, model_id: Optional[str] = None) -> None:
        """Сохраняет загруженную модель и сбрасывает её dirty-флаг."""
        mid = model_id or self.active_model_id
        model = self.models.get(mid) if mid else None
        if model is None or not model.loaded:
            return
        variant = model.variants.get("base")
        if variant is None or variant.composition is None:
            return
        assert mid is not None
        comp_node = variant.nodes.get("composition")
        correlations = (comp_node.params.get("correlations", {})
                        if comp_node is not None else {})
        self._repo.save_composition(mid, variant.composition, correlations)
        model.dirty = False
        self._notify(StateChange(StateChangeKind.MODEL_LIST))

    def save_all_dirty(self) -> list[str]:
        """Сохраняет все загруженные изменённые модели; возвращает их id."""
        saved: list[str] = []
        for mid, model in list(self.models.items()):
            if model.loaded and model.dirty:
                self.save_model(mid)
                saved.append(mid)
        return saved

    # --- узлы «Флэш» (история расчётов) ----------------------------------
    #
    # Расчёт выполняется во View в отдельном потоке (движок непрерываемый),
    # эти методы лишь заводят узел, переводят статус и хранят результат.

    def node_by_id(self, node_id: str) -> Optional[GraphNode]:
        variant = self.active_variant
        return None if variant is None else variant.nodes.get(node_id)

    def node_ref(self, node_id: str) -> Optional[NodeRef]:
        """Возвращает стабильный адрес узла активного варианта."""
        if (self.active_model_id is None or self.active_variant_id is None
                or self.node_by_id(node_id) is None):
            return None
        return NodeRef(self.active_model_id, self.active_variant_id, node_id)

    def node_by_ref(self, ref: NodeRef) -> Optional[GraphNode]:
        """Находит узел по полному адресу, не меняя активную модель."""
        model = self.models.get(ref.model_id)
        variant = model.variants.get(ref.variant_id) if model else None
        return variant.nodes.get(ref.node_id) if variant else None

    def composition_by_ref(self, ref: NodeRef) -> Optional[Composition]:
        """Возвращает состав варианта, которому принадлежит ``ref``."""
        model = self.models.get(ref.model_id)
        variant = model.variants.get(ref.variant_id) if model else None
        return variant.composition if variant else None

    def _resolve_node(self, node: str | NodeRef) -> Optional[GraphNode]:
        return self.node_by_ref(node) if isinstance(node, NodeRef) else self.node_by_id(node)

    def new_flash_run(self, p: Optional[float] = None,
                      t: Optional[float] = None) -> Optional[str]:
        """
        Заводит новый узел флэша (по умолчанию P/T как у последнего расчёта или
        пластовые) и открывает его вкладку. Возвращает id узла.
        """
        variant = self.active_variant
        composition = self.active_composition
        if variant is None or composition is None:
            return None
        self._push_undo()
        if p is None or t is None:
            runs = variant.flash_runs()
            if runs:
                p = runs[-1].params.get("P", 100.0) if p is None else p
                t = runs[-1].params.get("T", 20.0) if t is None else t
            else:
                p = 100.0 if p is None else p
                t = round(composition.T - 273.15, 2) if t is None else t

        variant.flash_seq += 1
        node_id = f"flash_{variant.flash_seq}"
        variant.nodes[node_id] = GraphNode(
            node_id=node_id,
            kind=NodeKind.FLASH,
            title="Flash",
            status=NodeStatus.EMPTY,
            params={"P": float(p), "T": float(t)},
            upstream=["composition"],
        )
        self.open_node(node_id)
        return node_id

    def new_experiment(self, kind: str, defaults: dict) -> Optional[str]:
        """
        Заводит новый узел-эксперимент (`kind` = cce/dle/separator) с
        параметрами по умолчанию и открывает его вкладку. Возвращает id узла.
        """
        variant = self.active_variant
        if variant is None or self.active_composition is None:
            return None
        self._push_undo()
        variant.exp_seq += 1
        node_id = f"exp_{variant.exp_seq}"
        params = {"kind": kind}
        params.update(defaults)
        variant.nodes[node_id] = GraphNode(
            node_id=node_id,
            kind=NodeKind.EXPERIMENT,
            title=kind.upper(),
            status=NodeStatus.EMPTY,
            params=params,
            upstream=["composition"],
        )
        self.open_node(node_id)
        return node_id

    def new_envelope(self, defaults: dict) -> Optional[str]:
        """
        Заводит новый узел фазовой огибающей (P-T + крит. точка + пластовое
        Psat) с параметрами по умолчанию и открывает его вкладку.
        Возвращает id узла.
        """
        variant = self.active_variant
        if variant is None or self.active_composition is None:
            return None
        self._push_undo()
        variant.env_seq += 1
        node_id = f"env_{variant.env_seq}"
        variant.nodes[node_id] = GraphNode(
            node_id=node_id,
            kind=NodeKind.PHASE_ENVELOPE,
            title="Phase envelope",
            status=NodeStatus.EMPTY,
            params=dict(defaults),
            upstream=["composition"],
        )
        self.open_node(node_id)
        return node_id

    def open_compare(self, member_ids: list[str]) -> Optional[str]:
        """
        Открывает вкладку сравнения выбранных узлов (единственный узел
        `compare`, его список участников обновляется). Нужно ≥2 существующих
        узла. Не откатывается undo (это представление, не правка данных).
        """
        variant = self.active_variant
        if variant is None:
            return None
        # упорядочиваем по порядку создания узлов (не по порядку выбора)
        wanted = set(member_ids)
        members = [nid for nid in variant.nodes if nid in wanted]
        if len(members) < 2:
            return None
        nid = "compare"
        if nid not in variant.nodes:
            variant.nodes[nid] = GraphNode(nid, NodeKind.COMPARE, "Compare",
                                           NodeStatus.FRESH, params={"members": members})
        else:
            variant.nodes[nid].params["members"] = members
        self.open_node(nid)
        return nid

    def duplicate_flash(self, node_id: str) -> Optional[str]:
        """Создаёт новый флэш-узел с теми же P/T (не запуская расчёт)."""
        node = self.node_by_id(node_id)
        if node is None or node.kind is not NodeKind.FLASH:
            return None
        return self.new_flash_run(node.params.get("P"), node.params.get("T"))

    def rename_node(self, node_id: str, name: str) -> None:
        """Задаёт (или сбрасывает пустым) пользовательскую подпись узла."""
        node = self.node_by_id(node_id)
        if node is None:
            return
        self._push_undo()
        name = (name or "").strip()
        if name:
            node.params["label"] = name
        else:
            node.params.pop("label", None)
        self._notify(StateChange(StateChangeKind.NODE, self.node_ref(node_id)))

    def delete_node(self, node_id: str) -> None:
        """Удаляет узел (кроме Composition) целиком — из истории и из вкладок."""
        variant = self.active_variant
        if variant is None:
            return
        node = variant.nodes.get(node_id)
        if node is None or node.kind is NodeKind.COMPOSITION:
            return
        self._push_undo()
        if node_id in variant.open_node_ids:
            idx = variant.open_node_ids.index(node_id)
            variant.open_node_ids.remove(node_id)
            if variant.active_node_id == node_id:
                variant.active_node_id = (
                    variant.open_node_ids[min(idx, len(variant.open_node_ids) - 1)]
                    if variant.open_node_ids else None)
        variant.nodes.pop(node_id, None)
        self._notify(StateChange(StateChangeKind.WORKSPACE))

    def _invalidate_flash(self) -> None:
        """Помечает посчитанные флэши/эксперименты/огибающие `STALE` при изменении состава."""
        variant = self.active_variant
        if variant is None:
            return
        for node in (variant.flash_runs() + variant.experiment_runs()
                     + variant.envelope_runs()):
            if node.status is NodeStatus.FRESH:
                node.status = NodeStatus.STALE

    def set_flash_params(self, node_id: str, p: float, t: float) -> None:
        """Сохраняет P (бар)/T (°C) в параметрах узла (без перерисовки)."""
        node = self.node_by_id(node_id)
        if node is not None:
            node.params["P"] = float(p)
            node.params["T"] = float(t)

    def update_node_params(self, node_id: str | NodeRef, params: dict,
                           *, notify: bool = True) -> None:
        """Обновляет параметры узла через state-инварианты, а не из View."""
        node = self._resolve_node(node_id)
        if node is None:
            return
        self._push_undo()
        node.params.update(params)
        if node.result is not None and node.status is NodeStatus.FRESH:
            node.status = NodeStatus.STALE
        if notify:
            ref = node_id if isinstance(node_id, NodeRef) else self.node_ref(node.node_id)
            self._notify(StateChange(StateChangeKind.NODE, ref))

    # --- фактические точки эксперимента ----------------------------------

    def _lab_node(self, node_id: str) -> Optional[GraphNode]:
        node = self.node_by_id(node_id)
        return node if node is not None and node.kind is NodeKind.EXPERIMENT else None

    @staticmethod
    def _ensure_lab_rows(
        node: GraphNode,
        columns: list[str],
    ) -> list[list[float | None]]:
        """Returns a rectangular, schema-compatible measured-data table."""
        data = node.params.get("lab_data")
        if not isinstance(data, dict) or data.get("columns") != columns:
            data = {"schema_version": 1, "columns": list(columns), "rows": []}
            node.params["lab_data"] = data
        raw_rows = data.get("rows")
        rows: list[list[float | None]] = []
        if isinstance(raw_rows, list):
            for raw_row in raw_rows:
                if isinstance(raw_row, list):
                    rows.append((raw_row[:len(columns)] + [None] * len(columns))
                                [:len(columns)])
                else:
                    rows.append([None] * len(columns))
        data["rows"] = rows
        return rows

    def add_lab_data_row(self, node_id: str, columns: list[str]) -> None:
        """Добавляет пустую строку фактических данных с поддержкой undo."""
        node = self._lab_node(node_id)
        if node is None or not columns:
            return
        self._push_undo()
        rows = self._ensure_lab_rows(node, columns)
        rows.append([None] * len(columns))
        self._mark_dirty()

    def append_lab_data_rows(self, node_id: str, columns: list[str],
                             new_rows: list[list[float | None]]) -> None:
        """Добавляет несколько валидированных строк фактических данных."""
        node = self._lab_node(node_id)
        if node is None or not columns or not new_rows:
            return
        self._push_undo()
        rows = self._ensure_lab_rows(node, columns)
        rows.extend([
            (row[:len(columns)] + [None] * len(columns))[:len(columns)]
            for row in new_rows
        ])
        self._mark_dirty()

    def paste_lab_data_rows(
        self,
        node_id: str,
        columns: list[str],
        start_row: int,
        new_rows: list[list[float | None]],
    ) -> None:
        """Writes pasted points from a cell, extending the table as needed.

        ``None`` values are left untouched so a partial clipboard table does
        not erase measurements in other columns.  This makes a one-column
        Excel paste land in the focused column while preserving the rest of
        each existing point.
        """
        node = self._lab_node(node_id)
        if node is None or not columns or not new_rows:
            return
        self._push_undo()
        rows = self._ensure_lab_rows(node, columns)
        start_row = max(0, int(start_row))
        while len(rows) < start_row + len(new_rows):
            rows.append([None] * len(columns))
        for row_offset, source_row in enumerate(new_rows):
            target_row = rows[start_row + row_offset]
            if not isinstance(target_row, list):
                target_row = [None] * len(columns)
                rows[start_row + row_offset] = target_row
            normalized = (source_row[:len(columns)]
                          + [None] * len(columns))[:len(columns)]
            for column_index, value in enumerate(normalized):
                if value is not None:
                    target_row[column_index] = value
        self._mark_dirty()

    def remove_lab_data_row(self, node_id: str) -> None:
        """Удаляет последнюю строку фактических данных."""
        node = self._lab_node(node_id)
        if node is None:
            return
        data = node.params.get("lab_data")
        columns = data.get("columns") if isinstance(data, dict) else None
        if (not isinstance(columns, list)
                or not all(isinstance(column, str) for column in columns)):
            return
        rows = self._ensure_lab_rows(node, columns)
        if rows:
            self._push_undo()
            rows.pop()
            self._mark_dirty()

    def clear_lab_data(self, node_id: str) -> None:
        """Очищает точки, сохраняя выбранную схему колонок."""
        node = self._lab_node(node_id)
        if node is None:
            return
        data = node.params.get("lab_data")
        if not isinstance(data, dict) or not isinstance(data.get("rows"), list):
            return
        if data["rows"]:
            self._push_undo()
            data["rows"] = []
            self._mark_dirty()

    def set_lab_data_value(self, node_id: str, row: int, column: int,
                           value: float | None) -> None:
        """Меняет одну фактическую точку без уведомления (сохраняет фокус поля)."""
        node = self._lab_node(node_id)
        if node is None:
            return
        data = node.params.get("lab_data")
        if not isinstance(data, dict):
            return
        columns = data.get("columns")
        if not isinstance(columns, list) or not columns:
            return
        rows = self._ensure_lab_rows(node, list(columns))
        if (row < 0 or row >= len(rows) or column < 0
                or column >= len(rows[row])):
            return
        if rows[row][column] == value:
            return
        self._push_undo()
        rows[row][column] = value
        self._mark_dirty()

    # --- обобщённые переходы статуса узла-расчёта (флэш/эксперимент) -----

    def set_node_running(self, node_id: str | NodeRef) -> None:
        """Переводит узел-расчёт в `RUNNING` и перерисовывает."""
        node = self._resolve_node(node_id)
        if node is not None:
            ref = node_id if isinstance(node_id, NodeRef) else self.node_ref(node.node_id)
            node.status = NodeStatus.RUNNING
            node.error = None
            self._notify(StateChange(StateChangeKind.NODE, ref))

    def set_node_result(self, node_id: str | NodeRef, result: object) -> None:
        """Сохраняет результат, узел → `FRESH`, перерисовка."""
        node = self._resolve_node(node_id)
        if node is not None:
            ref = node_id if isinstance(node_id, NodeRef) else self.node_ref(node.node_id)
            node.result = result
            node.status = NodeStatus.FRESH
            node.error = None
            self._notify(StateChange(StateChangeKind.NODE, ref))

    def set_node_error(self, node_id: str | NodeRef, message: str) -> None:
        """Сохраняет ошибку, узел → `STALE`, перерисовка."""
        node = self._resolve_node(node_id)
        if node is not None:
            ref = node_id if isinstance(node_id, NodeRef) else self.node_ref(node.node_id)
            node.result = None
            node.status = NodeStatus.STALE
            node.error = message
            self._notify(StateChange(StateChangeKind.NODE, ref))

    def reset_node(self, node_id: str | NodeRef) -> None:
        """Сбрасывает узел в `EMPTY` (например, при отмене), перерисовка."""
        node = self._resolve_node(node_id)
        if node is not None:
            ref = node_id if isinstance(node_id, NodeRef) else self.node_ref(node.node_id)
            node.status = NodeStatus.EMPTY
            node.error = None
            self._notify(StateChange(StateChangeKind.NODE, ref))

    # обратная совместимость (флэш-именованные вызовы)
    set_flash_running = set_node_running
    set_flash_result = set_node_result
    set_flash_error = set_node_error
    reset_flash = reset_node

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

    @property
    def active_node(self) -> Optional[GraphNode]:
        variant = self.active_variant
        if variant is None or variant.active_node_id is None:
            return None
        return variant.nodes.get(variant.active_node_id)
