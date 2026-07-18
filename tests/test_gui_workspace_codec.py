"""Workspace v3 сохраняет точные ids, подписи, stale-результаты и compare."""

from gui.app_state import GraphNode, NodeKind, NodeStatus, Variant
from gui.workspace_codec import restore_workspace, snapshot_workspace


def test_workspace_v3_roundtrip_preserves_nodes():
    variant = Variant("base", "Base")
    variant.nodes["composition"] = GraphNode(
        "composition", NodeKind.COMPOSITION, "Composition", NodeStatus.FRESH)
    variant.nodes["exp_7"] = GraphNode(
        "exp_7", NodeKind.EXPERIMENT, "DLE", NodeStatus.STALE,
        params={"kind": "dle", "label": "lab run"}, result={"rows": [[1.0]]},
        upstream=["composition"])
    variant.nodes["compare"] = GraphNode(
        "compare", NodeKind.COMPARE, "Compare", NodeStatus.FRESH,
        params={"members": ["exp_7", "other"]})
    variant.exp_seq = 7
    variant.open_node_ids = ["composition", "exp_7", "compare"]
    variant.active_node_id = "compare"

    snap = snapshot_workspace(variant)
    restored = Variant("base", "Base")
    restored.nodes["composition"] = variant.nodes["composition"]
    restore_workspace(restored, snap)

    assert restored.exp_seq == 7
    assert restored.nodes["exp_7"].params["label"] == "lab run"
    assert restored.nodes["exp_7"].status is NodeStatus.STALE
    assert restored.nodes["exp_7"].result == {"rows": [[1.0]]}
    assert restored.nodes["compare"].params["members"] == ["exp_7", "other"]
    assert restored.active_node_id == "compare"


def test_workspace_v2_is_migrated_in_memory():
    ws = {"flashes": [{"P": 50, "T": 20, "label": "old"}],
          "open_tabs": ["flash_1"], "active_tab": "flash_1"}
    variant = Variant("base", "Base")
    variant.nodes["composition"] = GraphNode(
        "composition", NodeKind.COMPOSITION, "Composition", NodeStatus.FRESH)
    restore_workspace(variant, ws)
    assert variant.flash_seq == 1
    assert variant.nodes["flash_1"].params["label"] == "old"
