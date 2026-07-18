"""Версия схемы, диагностика, backup/lock и concurrency ModelStore."""

import json

import pytest

from calc_core.Utils.Export import ModelJSONDB
from calc_core.Utils.ModelStore import (
    MODEL_SCHEMA_VERSION,
    ConcurrentModelStoreUpdate,
    ModelRecordValidationError,
    ModelStoreCorruptError,
    ModelStoreLockedError,
    model_store_lock,
    read_model_store,
    update_model_store,
    write_model_store,
)


def _record(name: str = "Fluid") -> dict:
    return {
        "Model_name": name,
        "composition": {"C1": 1.0},
        "composition_data": {},
        "eos": "PREOS",
        "results": [],
    }


def test_legacy_record_migrates_in_memory_and_on_next_write(tmp_path):
    path = tmp_path / "models.json"
    path.write_text(json.dumps({"A": _record()}), encoding="utf-8")

    data = read_model_store(path)
    assert data["A"]["schema_version"] == MODEL_SCHEMA_VERSION
    assert "schema_version" not in path.read_text(encoding="utf-8")

    write_model_store(path, data)
    saved = json.loads(path.read_text(encoding="utf-8"))
    assert saved["A"]["schema_version"] == MODEL_SCHEMA_VERSION
    assert (tmp_path / "models.json.bak").exists()


def test_corrupt_json_reports_location(tmp_path):
    path = tmp_path / "models.json"
    path.write_text('{"A":', encoding="utf-8")
    with pytest.raises(ModelStoreCorruptError, match="строка 1"):
        read_model_store(path)


@pytest.mark.parametrize(
    "record, message",
    [
        ({}, "composition"),
        ({**_record(), "schema_version": 999}, "schema_version=999"),
        ({**_record(), "T_res": float("nan")}, "T_res"),
        ({**_record(), "results": {}}, "results"),
        ({**_record(), "composition_data": {"Pc": {"C1": float("nan")}}},
         "сериализовать"),
        ({**_record(), "Model_name": 42}, "Model_name"),
    ],
)
def test_invalid_record_has_model_context(record, message, tmp_path):
    path = tmp_path / "models.json"
    path.write_text(json.dumps({"BROKEN": record}), encoding="utf-8")
    with pytest.raises(ModelRecordValidationError, match=message) as exc:
        read_model_store(path)
    assert "BROKEN" in str(exc.value)


def test_backups_rotate_for_each_successful_update(tmp_path):
    path = tmp_path / "models.json"
    write_model_store(path, {"A": _record("v0")})

    for title in ("v1", "v2", "v3"):
        update_model_store(
            path,
            lambda data, value=title: data["A"].__setitem__("Model_name", value),
        )

    assert read_model_store(path)["A"]["Model_name"] == "v3"
    assert read_model_store(str(path) + ".bak")["A"]["Model_name"] == "v2"
    assert read_model_store(str(path) + ".bak.1")["A"]["Model_name"] == "v1"
    assert read_model_store(str(path) + ".bak.2")["A"]["Model_name"] == "v0"


def test_model_json_db_refuses_stale_snapshot(tmp_path):
    path = str(tmp_path / "models.json")
    seed = ModelJSONDB(path)
    seed.export("A", "A", {"C1": 1.0}, {}, "PREOS")
    seed.save()

    first = ModelJSONDB(path)
    stale = ModelJSONDB(path)
    first.export("B", "B", {"C1": 1.0}, {}, "PREOS")
    first.save()
    stale.export("C", "C", {"C1": 1.0}, {}, "PREOS")
    with pytest.raises(ConcurrentModelStoreUpdate):
        stale.save()
    assert set(read_model_store(path)) == {"A", "B"}


def test_live_lock_is_not_silently_broken(tmp_path):
    path = tmp_path / "models.json"
    lock = tmp_path / "models.json.lock"
    lock.write_text("pid=test\n", encoding="ascii")
    with pytest.raises(ModelStoreLockedError):
        with model_store_lock(path, timeout=0.01, stale_after=3600):
            pass
    assert lock.exists()


def test_lock_cleanup_does_not_remove_replacement_owner(tmp_path):
    path = tmp_path / "models.json"
    lock = tmp_path / "models.json.lock"
    with model_store_lock(path):
        lock.write_text("pid=replacement\n", encoding="ascii")
    assert lock.read_text(encoding="ascii") == "pid=replacement\n"
