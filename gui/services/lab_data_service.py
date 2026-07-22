"""Project- and model-scoped laboratory data catalog.

``models.json`` remains a strict ``model_id -> model`` mapping consumed by the
calculation layer.  The GUI catalog deliberately lives in the adjacent
``models.lab.json`` file, while experiment nodes keep only ``lab_data_ref``.
This avoids duplicating measured points across model variants and keeps the
calculation database backward-compatible.
"""

from __future__ import annotations

import json
import math
import shutil
from datetime import datetime
from pathlib import Path
from typing import Any
from uuid import uuid4

from calc_core.Utils.AtomicFile import atomic_write_json
from calc_core.Utils.ModelStore import model_store_lock

_SCHEMA_VERSION = 1
_SCOPES = {"project", "model"}


class LabDataStoreError(RuntimeError):
    """The auxiliary Lab Data catalog is unreadable or incompatible."""


def catalog_path(models_db_path: str | Path) -> Path:
    """Returns the GUI-only project catalog stored next to ``models.json``."""
    return Path(models_db_path).with_suffix(".lab.json")


def _empty_store() -> dict[str, Any]:
    return {"schema_version": _SCHEMA_VERSION, "projects": {}}


def _finite(value: object) -> float | None:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return None
    number = float(value)
    return number if math.isfinite(number) else None


def _normalise_dataset(raw: object) -> dict[str, Any] | None:
    if not isinstance(raw, dict):
        return None
    dataset_id = raw.get("dataset_id")
    kind = raw.get("experiment_kind")
    columns = raw.get("columns")
    scope = raw.get("scope", "project")
    if (not isinstance(dataset_id, str) or not dataset_id
            or not isinstance(kind, str) or not kind
            or not isinstance(columns, list)
            or not all(isinstance(column, str) and column for column in columns)
            or not isinstance(scope, str) or scope not in _SCOPES):
        return None
    model_id = raw.get("model_id")
    if scope == "model" and (not isinstance(model_id, str) or not model_id):
        return None
    rows: list[list[float | None]] = []
    for raw_row in raw.get("rows", []):
        if not isinstance(raw_row, list):
            continue
        rows.append([
            _finite(value) for value in
            (raw_row[:len(columns)] + [None] * len(columns))[:len(columns)]
        ])
    raw_conditions = raw.get("conditions")
    conditions = dict(raw_conditions) if isinstance(raw_conditions, dict) else {}
    return {
        "dataset_id": dataset_id,
        "title": str(raw.get("title") or kind.upper()),
        "experiment_kind": kind.lower(),
        "scope": scope,
        **({"model_id": model_id} if scope == "model" else {}),
        "columns": list(columns),
        "rows": rows,
        "conditions": conditions,
        "created_at": str(raw.get("created_at") or ""),
    }


def _read(path: Path) -> dict[str, Any]:
    if not path.exists():
        return _empty_store()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise LabDataStoreError(f"Could not read Lab Data catalog: {exc}") from exc
    if not isinstance(raw, dict) or not isinstance(raw.get("projects"), dict):
        raise LabDataStoreError("Lab Data catalog has an invalid root structure")
    projects: dict[str, dict[str, list[dict[str, Any]]]] = {}
    for project_id, project in raw["projects"].items():
        if not isinstance(project_id, str) or not isinstance(project, dict):
            raise LabDataStoreError("Lab Data catalog has an invalid project entry")
        raw_datasets = project.get("datasets", [])
        if not isinstance(raw_datasets, list):
            raise LabDataStoreError("Lab Data catalog has an invalid dataset list")
        datasets = []
        for raw_dataset in raw_datasets:
            dataset = _normalise_dataset(raw_dataset)
            if dataset is None:
                raise LabDataStoreError("Lab Data catalog has an invalid dataset entry")
            datasets.append(dataset)
        projects[project_id] = {"datasets": datasets}
    return {"schema_version": _SCHEMA_VERSION, "projects": projects}


def _update(models_db_path: str | Path, update) -> None:
    path = catalog_path(models_db_path)
    with model_store_lock(path):
        data = _read(path)
        update(data)
        if path.exists():
            shutil.copy2(path, Path(str(path) + ".bak"))
        atomic_write_json(path, data)


def list_datasets(
    models_db_path: str | Path,
    project_id: str,
    *,
    experiment_kind: str | None = None,
    model_id: str | None = None,
) -> list[dict[str, Any]]:
    """Lists shared datasets and, when requested, datasets local to a model."""
    project = _read(catalog_path(models_db_path))["projects"].get(project_id, {})
    result = []
    for dataset in project.get("datasets", []):
        if experiment_kind and dataset["experiment_kind"] != experiment_kind.lower():
            continue
        if dataset["scope"] == "model" and dataset.get("model_id") != model_id:
            continue
        result.append(dict(dataset))
    return result


def get_dataset(models_db_path: str | Path, project_id: str,
                dataset_id: str, *, model_id: str | None = None) -> dict[str, Any] | None:
    for dataset in list_datasets(models_db_path, project_id, model_id=model_id):
        if dataset["dataset_id"] == dataset_id:
            return dataset
    return None


def create_dataset(
    models_db_path: str | Path,
    project_id: str,
    *,
    title: str,
    experiment_kind: str,
    columns: list[str],
    rows: list[list[float | None]],
    scope: str = "project",
    model_id: str | None = None,
    conditions: dict[str, float] | None = None,
) -> dict[str, Any]:
    """Atomically publishes a local table to the shared/model catalog."""
    scope = scope if scope in _SCOPES else "project"
    if scope == "model" and not model_id:
        raise ValueError("model-scoped Lab Data requires model_id")
    dataset = _normalise_dataset({
        "dataset_id": f"lab_{uuid4().hex}",
        "title": title.strip() or experiment_kind.upper(),
        "experiment_kind": experiment_kind,
        "scope": scope,
        "model_id": model_id,
        "columns": columns,
        "rows": rows,
        "conditions": conditions or {},
        "created_at": datetime.now().isoformat(timespec="seconds"),
    })
    if dataset is None:
        raise ValueError("invalid Lab Data dataset")

    def append(data) -> None:
        projects = data.setdefault("projects", {})
        project = projects.setdefault(project_id, {"datasets": []})
        project.setdefault("datasets", []).append(dataset)

    _update(models_db_path, append)
    return dict(dataset)


def update_dataset(
    models_db_path: str | Path,
    project_id: str,
    dataset_id: str,
    *,
    title: str,
    columns: list[str],
    rows: list[list[float | None]],
    conditions: dict[str, float] | None = None,
    model_id: str | None = None,
) -> dict[str, Any] | None:
    """Updates only an existing dataset visible in the requested scope."""
    updated: dict[str, Any] | None = None

    def replace(data) -> None:
        nonlocal updated
        project = data.setdefault("projects", {}).get(project_id)
        if not isinstance(project, dict):
            return
        for index, raw in enumerate(project.get("datasets", [])):
            dataset = _normalise_dataset(raw)
            if (dataset is None or dataset["dataset_id"] != dataset_id
                    or (dataset["scope"] == "model"
                        and dataset.get("model_id") != model_id)):
                continue
            candidate = _normalise_dataset({
                **dataset, "title": title, "columns": columns, "rows": rows,
                "conditions": conditions or {},
            })
            if candidate is None:
                raise ValueError("invalid Lab Data dataset")
            project["datasets"][index] = candidate
            updated = dict(candidate)
            return

    _update(models_db_path, replace)
    return updated


def delete_dataset(models_db_path: str | Path, project_id: str,
                   dataset_id: str, *, model_id: str | None = None) -> bool:
    deleted = False

    def remove(data) -> None:
        nonlocal deleted
        project = data.setdefault("projects", {}).get(project_id)
        if not isinstance(project, dict):
            return
        kept = []
        for raw in project.get("datasets", []):
            dataset = _normalise_dataset(raw)
            if (dataset is not None and dataset["dataset_id"] == dataset_id
                    and (dataset["scope"] != "model"
                         or dataset.get("model_id") == model_id)):
                deleted = True
                continue
            kept.append(raw)
        project["datasets"] = kept

    _update(models_db_path, remove)
    return deleted


def delete_project_datasets(models_db_path: str | Path, project_id: str) -> None:
    if not catalog_path(models_db_path).exists():
        return

    def remove(data) -> None:
        data.setdefault("projects", {}).pop(project_id, None)
    _update(models_db_path, remove)


def delete_model_datasets(models_db_path: str | Path, project_id: str,
                          model_id: str) -> None:
    if not catalog_path(models_db_path).exists():
        return

    def remove(data) -> None:
        project = data.setdefault("projects", {}).get(project_id)
        if isinstance(project, dict):
            project["datasets"] = [
                dataset for dataset in project.get("datasets", [])
                if not (dataset.get("scope") == "model"
                        and dataset.get("model_id") == model_id)
            ]
    _update(models_db_path, remove)


def effective_lab_data(models_db_path: str | Path, project_id: str,
                       model_id: str, params: dict) -> dict[str, Any] | None:
    """Resolves an experiment reference, falling back to its legacy local table."""
    ref = params.get("lab_data_ref") if isinstance(params, dict) else None
    if isinstance(ref, str) and ref:
        try:
            dataset = get_dataset(models_db_path, project_id, ref, model_id=model_id)
        except LabDataStoreError:
            dataset = None
        if dataset is not None:
            return {"columns": dataset["columns"], "rows": dataset["rows"],
                    "title": dataset["title"], "scope": dataset["scope"]}
    data = params.get("lab_data") if isinstance(params, dict) else None
    return data if isinstance(data, dict) else None
