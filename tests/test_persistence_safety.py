"""Атомарная запись и исправленные механические дефекты persist-слоя."""

from calc_core.Utils.AtomicFile import atomic_write_json
from calc_core.Utils.Results import ResultStore


def test_atomic_json_replaces_target(tmp_path):
    path = tmp_path / "data.json"
    path.write_text('{"old": true}', encoding="utf-8")
    atomic_write_json(path, {"new": 1})
    assert path.read_text(encoding="utf-8").strip().startswith('{')
    assert '"new": 1' in path.read_text(encoding="utf-8")
    assert not list(tmp_path.glob("*.tmp"))


def test_result_store_save_writes_pickle(tmp_path):
    store = ResultStore()
    store.add("test", {"x": 1}, {"y": 2})
    path = tmp_path / "result.pkl"
    store.save(str(path))
    assert path.stat().st_size > 0
