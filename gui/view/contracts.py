"""Явные контракты между composition root и feature-view.

Feature-view пока подключаются через mixin-композицию, но получают состояние,
сессию и координатор задач только через один типизированный ``ViewContext``.
Это переходный слой к объектной композиции без скрытых ``Any``-зависимостей.
"""

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from gui.app_state import AppState, GraphNode
from gui.calculation_coordinator import CalculationCoordinator
from gui.session import SessionState

DpgId = int | str


@runtime_checkable
class ViewHost(Protocol):
    """Минимальный общий API, который composition root даёт feature-view."""

    def _set_status(self, text: str) -> None: ...

    def _theme_stale(self) -> DpgId: ...

    def _track_modal(self, win: int) -> int: ...

    def _close_tracked_modal(self, win: int) -> None: ...


@dataclass(frozen=True, slots=True)
class ViewContext:
    """Долгоживущие framework-independent зависимости всех view-срезов."""

    state: AppState
    session: SessionState
    jobs: CalculationCoordinator


class ContextBoundView:
    """Типизированная база feature-view, привязанных к ``ViewContext``."""

    _view_context: ViewContext

    @property
    def _state(self) -> AppState:
        return self._view_context.state

    @property
    def _session(self) -> SessionState:
        return self._view_context.session

    @property
    def _jobs(self) -> CalculationCoordinator:
        return self._view_context.jobs

    # Методы ниже реализует PVTcalcApp либо один из feature-view. Заглушки
    # задают проверяемый контракт и находятся последними в MRO.
    def _set_status(self, text: str) -> None:
        raise NotImplementedError

    def _theme_stale(self) -> DpgId:
        raise NotImplementedError

    def _track_modal(self, win: int) -> int:
        raise NotImplementedError

    def _close_tracked_modal(self, win: int) -> None:
        raise NotImplementedError

    def _copy_table(self, columns, rows, label: str = "Table") -> None:
        raise NotImplementedError

    def _fmt(self, value: object) -> str:
        raise NotImplementedError

    def _g(self, value: object) -> str:
        raise NotImplementedError

    def _format_saved_at(self, iso: str | None) -> str:
        raise NotImplementedError

    def _restore_workspace(self, model_id: str) -> None:
        raise NotImplementedError

    def _arm_flash_poll(self) -> None:
        raise NotImplementedError

    def _schedule_session_autosave(self) -> None:
        raise NotImplementedError

    def _on_flash_cancel(self, sender, app_data, user_data) -> None:
        raise NotImplementedError

    def _on_duplicate_model_confirm(self, sender, app_data, user_data) -> None:
        raise NotImplementedError

    def _on_delete_model_confirm(self, sender, app_data, user_data) -> None:
        raise NotImplementedError

    def _on_view_composition(self, sender, app_data, user_data) -> None:
        raise NotImplementedError

    def _on_open_model_report_dialog(self, sender=None, app_data=None,
                                     user_data=None) -> None:
        raise NotImplementedError

    def _open_envelope_dialog(self, node_id: str | None) -> None:
        raise NotImplementedError

    def _ctrl_down(self) -> bool:
        raise NotImplementedError

    def _render_tree(self) -> None:
        raise NotImplementedError

    def _render_workspace(self) -> None:
        raise NotImplementedError

    def _open_lab_dataset_editor(self, dataset) -> None:
        raise NotImplementedError

    def _overview_is_selected(self) -> bool:
        raise NotImplementedError

    def _chart_grid_columns(self) -> int:
        raise NotImplementedError

    def _chart_card_width(self, columns: int) -> int:
        raise NotImplementedError

    def _chart_card_height(self) -> int:
        raise NotImplementedError

    def _chart_plot_height(self) -> int:
        raise NotImplementedError

    def _envelope_plot_height(self) -> int:
        raise NotImplementedError

    def _render_composition_tab(self, parent: DpgId, node: GraphNode) -> None:
        raise NotImplementedError

    def _render_flash_tab(self, parent: DpgId, node: GraphNode) -> None:
        raise NotImplementedError

    def _render_compare_tab(self, parent: DpgId, node: GraphNode) -> None:
        raise NotImplementedError

    def _render_experiment_tab(self, parent: DpgId, node: GraphNode) -> None:
        raise NotImplementedError

    def _render_envelope_tab(self, parent: DpgId, node: GraphNode) -> None:
        raise NotImplementedError
