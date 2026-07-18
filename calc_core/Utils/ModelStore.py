"""Версионированное и безопасное хранилище записей моделей.

Формат файла остаётся обратно совместимым: корень — словарь ``model_id ->
record``. Версия хранится внутри каждой записи, поэтому старый proxy
``Composition.from_db`` и внешние скрипты по-прежнему видят только модели.
"""

from __future__ import annotations

import json
import math
import os
import secrets
import shutil
import time
from collections.abc import Callable, Iterator, MutableMapping
from contextlib import contextmanager
from pathlib import Path
from typing import Any, TypeAlias

from calc_core.Utils.AtomicFile import atomic_write_json
from calc_core.Utils.Errors import PVTCalcError

MODEL_SCHEMA_VERSION = 1
DEFAULT_BACKUP_COUNT = 3

ModelRecord: TypeAlias = dict[str, Any]
ModelStoreData: TypeAlias = dict[str, ModelRecord]
StoreRevision: TypeAlias = tuple[int, int] | None


class ModelStoreError(PVTCalcError):
    """Базовая ожидаемая ошибка хранилища моделей."""


class ModelStoreCorruptError(ModelStoreError):
    """JSON повреждён либо корневая структура имеет неверный формат."""


class ModelRecordValidationError(ModelStoreError):
    """Одна из записей модели не соответствует поддерживаемой схеме."""


class ModelStoreLockedError(ModelStoreError):
    """Другой процесс не освободил lock-файл за допустимое время."""


class ConcurrentModelStoreUpdate(ModelStoreError):
    """Файл изменился после чтения и не может быть безопасно перезаписан."""


def store_revision(path: str | Path) -> StoreRevision:
    """Возвращает дешёвый fingerprint файла для optimistic concurrency."""
    target = Path(path)
    if not target.exists():
        return None
    stat = target.stat()
    return stat.st_mtime_ns, stat.st_size


def _validation_error(model_id: str, message: str) -> ModelRecordValidationError:
    return ModelRecordValidationError(f"Модель '{model_id}': {message}")


def validate_model_record(model_id: str, record: object) -> ModelRecord:
    """Проверяет и в памяти мигрирует одну запись к текущей схеме."""
    if not isinstance(model_id, str) or not model_id or not model_id.isidentifier():
        raise ModelRecordValidationError(
            f"Некорректный model_id {model_id!r}: требуется Python identifier")
    if not isinstance(record, dict):
        raise _validation_error(model_id, "запись должна быть JSON-объектом")

    version = record.get("schema_version", MODEL_SCHEMA_VERSION)
    if isinstance(version, bool) or not isinstance(version, int):
        raise _validation_error(model_id, "schema_version должен быть целым числом")
    if version < 1 or version > MODEL_SCHEMA_VERSION:
        raise _validation_error(
            model_id,
            f"неподдерживаемая schema_version={version}; "
            f"поддерживается {MODEL_SCHEMA_VERSION}",
        )

    composition = record.get("composition")
    if not isinstance(composition, dict) or not composition:
        raise _validation_error(model_id, "composition должен быть непустым объектом")
    total = 0.0
    for component, value in composition.items():
        if not isinstance(component, str) or not component:
            raise _validation_error(model_id, "имя компонента должно быть строкой")
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise _validation_error(model_id, f"доля {component} должна быть числом")
        number = float(value)
        if not math.isfinite(number) or number <= 0:
            raise _validation_error(
                model_id, f"доля {component} должна быть конечной и > 0")
        total += number
    if not math.isfinite(total) or total <= 0:
        raise _validation_error(model_id, "сумма состава должна быть конечной и > 0")

    composition_data = record.get("composition_data")
    if not isinstance(composition_data, dict):
        raise _validation_error(model_id, "composition_data должен быть объектом")
    results = record.get("results", [])
    if not isinstance(results, list):
        raise _validation_error(model_id, "results должен быть массивом")
    correlations = record.get("correlations")
    if correlations is not None and not isinstance(correlations, dict):
        raise _validation_error(model_id, "correlations должен быть объектом")
    for key in ("Model_name", "Field", "eos", "project_id", "created_at",
                "updated_at"):
        value = record.get(key)
        if value is not None and not isinstance(value, str):
            raise _validation_error(model_id, f"{key} должен быть строкой")
    t_res = record.get("T_res")
    if t_res is not None:
        if isinstance(t_res, bool) or not isinstance(t_res, (int, float)):
            raise _validation_error(model_id, "T_res должен быть числом")
        if not math.isfinite(float(t_res)) or float(t_res) <= 0:
            raise _validation_error(model_id, "T_res должен быть конечным и > 0 K")

    # Миграция v0 -> v1 выполняется только в загруженном объекте. На диск она
    # попадёт при следующем реальном сохранении пользователем.
    record["schema_version"] = MODEL_SCHEMA_VERSION
    record.setdefault("results", [])
    try:
        # Стандартный JSON не поддерживает NaN/Infinity. Проверяем всю запись,
        # включая composition_data/results, до ротации backup и записи файла.
        json.dumps(record, allow_nan=False)
    except (TypeError, ValueError) as exc:
        raise _validation_error(
            model_id, f"запись нельзя сериализовать в JSON: {exc}") from exc
    return record


def validate_model_store(data: object) -> ModelStoreData:
    """Валидирует корень и все записи; возвращает тот же mutable store."""
    if not isinstance(data, dict):
        raise ModelStoreCorruptError("Корень models.json должен быть JSON-объектом")
    for model_id, record in data.items():
        validate_model_record(model_id, record)
    return data


def read_model_store(path: str | Path, *, missing_ok: bool = False) -> ModelStoreData:
    """Читает, диагностирует JSON и валидирует все записи моделей."""
    target = Path(path)
    if not target.exists():
        if missing_ok:
            return {}
        raise FileNotFoundError(f"Файл моделей не найден: {target}")
    try:
        with target.open("r", encoding="utf-8") as stream:
            raw = json.load(stream)
    except json.JSONDecodeError as exc:
        raise ModelStoreCorruptError(
            f"Повреждён JSON {target}: строка {exc.lineno}, "
            f"столбец {exc.colno}: {exc.msg}") from exc
    except OSError as exc:
        raise ModelStoreError(f"Не удалось прочитать {target}: {exc}") from exc
    return validate_model_store(raw)


def _lock_path(path: Path) -> Path:
    return Path(str(path) + ".lock")


@contextmanager
def model_store_lock(
    path: str | Path,
    *,
    timeout: float = 3.0,
    stale_after: float = 60.0,
) -> Iterator[None]:
    """Межпроцессный lock на базе атомарного создания соседнего файла."""
    target = Path(path)
    lock = _lock_path(target)
    owner_token = f"pid={os.getpid()};token={secrets.token_hex(12)}\n"
    deadline = time.monotonic() + timeout
    lock.parent.mkdir(parents=True, exist_ok=True)
    while True:
        try:
            descriptor = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
            try:
                os.write(descriptor, owner_token.encode("ascii"))
            finally:
                os.close(descriptor)
            break
        except FileExistsError:
            try:
                age = time.time() - lock.stat().st_mtime
                if age > stale_after:
                    lock.unlink(missing_ok=True)
                    continue
            except FileNotFoundError:
                continue
            if time.monotonic() >= deadline:
                raise ModelStoreLockedError(
                    f"Хранилище занято другим процессом: {lock}")
            time.sleep(0.05)
    try:
        yield
    finally:
        # Если наш lock признали протухшим и заменили, нельзя удалять lock
        # нового владельца при выходе исходного процесса.
        try:
            if lock.read_text(encoding="ascii") == owner_token:
                lock.unlink(missing_ok=True)
        except FileNotFoundError:
            pass


def _rotate_backups(path: Path, count: int) -> None:
    if count <= 0 or not path.exists():
        return
    for index in range(count - 1, 0, -1):
        source = Path(str(path) + (".bak" if index == 1 else f".bak.{index - 1}"))
        destination = Path(str(path) + f".bak.{index}")
        if source.exists():
            shutil.copy2(source, destination)
    shutil.copy2(path, Path(str(path) + ".bak"))


_NO_REVISION_CHECK = object()


def _write_unlocked(
    path: Path,
    data: ModelStoreData,
    *,
    backup_count: int,
    expected_revision: StoreRevision | object,
) -> None:
    if (expected_revision is not _NO_REVISION_CHECK
            and store_revision(path) != expected_revision):
        raise ConcurrentModelStoreUpdate(
            f"{path} изменён другим процессом после чтения; сохранение отменено")
    validate_model_store(data)
    _rotate_backups(path, backup_count)
    atomic_write_json(path, data)


def write_model_store(
    path: str | Path,
    data: ModelStoreData,
    *,
    backup_count: int = DEFAULT_BACKUP_COUNT,
    expected_revision: StoreRevision | object = _NO_REVISION_CHECK,
) -> None:
    """Под lock валидирует, ротирует backup и атомарно пишет весь store."""
    target = Path(path)
    with model_store_lock(target):
        _write_unlocked(
            target,
            data,
            backup_count=backup_count,
            expected_revision=expected_revision,
        )


def update_model_store(
    path: str | Path,
    update: Callable[[MutableMapping[str, ModelRecord]], None],
    *,
    missing_ok: bool = False,
    backup_count: int = DEFAULT_BACKUP_COUNT,
) -> ModelStoreData:
    """Атомарная read-modify-write транзакция под одним lock."""
    target = Path(path)
    with model_store_lock(target):
        data = read_model_store(target, missing_ok=missing_ok)
        update(data)
        _write_unlocked(
            target,
            data,
            backup_count=backup_count,
            expected_revision=_NO_REVISION_CHECK,
        )
        return data
