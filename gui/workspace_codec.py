"""Сериализация рабочего пространства независимо от DearPyGui."""

from __future__ import annotations

import logging
from typing import Any

from gui.app_state import GraphNode, NodeKind, NodeStatus, Variant
from gui.services import flash_service

logger = logging.getLogger(__name__)


def _result_snapshot(node: GraphNode):
    if node.result is None:
        return None
    if node.kind is NodeKind.FLASH:
        return flash_service.snapshot_flash_result(node.result)
    return node.result


def snapshot_workspace(variant: Variant) -> dict:
    """Возвращает полное JSON-совместимое состояние workspace v3."""
    nodes = []
    for node in variant.nodes.values():
        if node.kind is NodeKind.COMPOSITION:
            continue
        try:
            result = _result_snapshot(node)
        except Exception:  # noqa: BLE001
            logger.warning("Не удалось сериализовать результат узла %s", node.node_id)
            result = None
        nodes.append({
            "node_id": node.node_id,
            "kind": node.kind.name,
            "title": node.title,
            "status": node.status.name,
            "params": dict(node.params),
            "result": result,
            "error": node.error,
            "upstream": list(node.upstream),
        })
    return {
        "schema_version": 3,
        "open_tabs": list(variant.open_node_ids),
        "active_tab": variant.active_node_id,
        "sequences": {
            "flash": variant.flash_seq,
            "experiment": variant.exp_seq,
            "envelope": variant.env_seq,
        },
        "nodes": nodes,
    }


def _status(value: Any, *, has_result: bool) -> NodeStatus:
    try:
        status = NodeStatus[str(value)]
    except (KeyError, TypeError):
        status = NodeStatus.FRESH if has_result else NodeStatus.EMPTY
    # RUNNING нельзя восстановить после перезапуска процесса.
    return NodeStatus.STALE if status is NodeStatus.RUNNING else status


def _restore_result(kind: NodeKind, snap):
    if snap is None:
        return None
    return flash_service.restore_flash_result(snap) if kind is NodeKind.FLASH else snap


def _legacy_nodes(ws: dict) -> list[dict]:
    """Преобразует workspace v2 в общий список узлов v3."""
    out: list[dict] = []
    flashes = ws.get("flashes", [])
    if not isinstance(flashes, list):
        flashes = []
    for i, item in enumerate((x for x in flashes if isinstance(x, dict)), start=1):
        out.append({
            "node_id": f"flash_{i}", "kind": "FLASH", "title": "Flash",
            "params": {"P": item.get("P"), "T": item.get("T"),
                       **({"label": item["label"]} if item.get("label") else {})},
            "result": item.get("result"), "upstream": ["composition"],
        })
    experiments = ws.get("experiments", [])
    if not isinstance(experiments, list):
        experiments = []
    for i, item in enumerate((x for x in experiments if isinstance(x, dict)), start=1):
        params = dict(item.get("params", {}))
        params.setdefault("kind", item.get("kind"))
        out.append({
            "node_id": f"exp_{i}", "kind": "EXPERIMENT",
            "title": str(item.get("kind") or "Experiment").upper(),
            "params": params, "result": item.get("result"),
            "upstream": ["composition"],
        })
    envelopes = ws.get("envelopes", [])
    if not isinstance(envelopes, list):
        envelopes = []
    for i, item in enumerate((x for x in envelopes if isinstance(x, dict)), start=1):
        out.append({
            "node_id": f"env_{i}", "kind": "PHASE_ENVELOPE",
            "title": "Phase envelope", "params": dict(item.get("params", {})),
            "result": item.get("result"), "upstream": ["composition"],
        })
    return out


def restore_workspace(variant: Variant, ws: dict) -> None:
    """Восстанавливает workspace v2/v3 в уже загруженный ``variant``."""
    if not isinstance(ws, dict):
        logger.warning("Пропущен workspace с неверной структурой: %r", type(ws).__name__)
        return
    raw_records = ws.get("nodes")
    records: list[dict] = raw_records if isinstance(raw_records, list) else _legacy_nodes(ws)
    for rec in records:
        if not isinstance(rec, dict):
            logger.warning("Пропущена запись узла неверного типа: %r", type(rec).__name__)
            continue
        try:
            kind = NodeKind[str(rec.get("kind"))]
        except (KeyError, TypeError):
            logger.warning("Неизвестный тип узла в сессии: %s", rec.get("kind"))
            continue
        if kind is NodeKind.COMPOSITION:
            continue
        node_id = str(rec.get("node_id") or "")
        if not node_id:
            continue
        try:
            result = _restore_result(kind, rec.get("result"))
        except Exception:  # noqa: BLE001
            logger.warning("Не удалось восстановить результат узла %s", node_id)
            result = None
        raw_params = rec.get("params")
        params = dict(raw_params) if isinstance(raw_params, dict) else {}
        raw_upstream = rec.get("upstream")
        upstream = [item for item in raw_upstream if isinstance(item, str)] \
            if isinstance(raw_upstream, list) else []
        error = rec.get("error")
        variant.nodes[node_id] = GraphNode(
            node_id=node_id,
            kind=kind,
            title=str(rec.get("title") or kind.name.title()),
            status=_status(rec.get("status"), has_result=result is not None),
            params=params,
            result=result,
            error=error if isinstance(error, str) else None,
            upstream=upstream,
        )

    seq = ws.get("sequences") or {}
    if not isinstance(seq, dict):
        seq = {}

    def sequence(value) -> int:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            return 0
        try:
            return max(0, int(value))
        except (ValueError, OverflowError):
            return 0

    variant.flash_seq = max(sequence(seq.get("flash")), _max_suffix(variant, "flash_"))
    variant.exp_seq = max(sequence(seq.get("experiment")), _max_suffix(variant, "exp_"))
    variant.env_seq = max(sequence(seq.get("envelope")), _max_suffix(variant, "env_"))
    open_tabs = ws.get("open_tabs")
    if isinstance(open_tabs, list):
        variant.open_node_ids = [nid for nid in open_tabs
                                 if isinstance(nid, str) and nid in variant.nodes]
    active = ws.get("active_tab")
    variant.active_node_id = (active if isinstance(active, str) and active in variant.nodes else
                              (variant.open_node_ids[-1] if variant.open_node_ids else None))


def _max_suffix(variant: Variant, prefix: str) -> int:
    values = []
    for node_id in variant.nodes:
        if node_id.startswith(prefix):
            try:
                values.append(int(node_id[len(prefix):]))
            except ValueError:
                pass
    return max(values, default=0)
