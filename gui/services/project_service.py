"""
Сервис страницы Projects: каталог компонентов БД, создание нового флюида,
импорт состава из Excel, сводка «что рассчитано» по модели.

Не импортирует DearPyGui. Опирается на движок:
- `JsonDBReader` — каталог доступных компонентов (`DB.json`);
- `Composition(zi, T_res, eos_name)` + `evaluate_composition_data` — создание;
- `ModelJSONDB.export_and_save(..., t_res=...)` — сохранение в models.json
  (T_res — новое опциональное поле записи, см. Export.py);
- `CompositionExcelLoader` — импорт zi из Excel (2 колонки: имя | доля).
"""

import logging
import math
import re
from copy import deepcopy
from datetime import datetime
from pathlib import Path
from typing import Optional

from calc_core.Composition.Composition import Composition
from calc_core.EOS.BaseEOS import EOSType
from calc_core.Utils.CompositionLoader import CompositionExcelLoader
from calc_core.Utils.E300Import import parse_e300
from calc_core.Utils.Export import ModelJSONDB
from calc_core.Utils.JsonDBReader import JsonDBReader
from calc_core.Utils.ModelStore import read_model_store, update_model_store
from gui.services import composition_service as comp_svc

logger = logging.getLogger(__name__)

# кэш БД компонентов (читается один раз за процесс)
_db_cache: Optional[dict] = None


def _component_db() -> dict:
    global _db_cache
    if _db_cache is None:
        _db_cache = JsonDBReader().load_database('DB.json')
    return _db_cache


# --- каталог компонентов -----------------------------------------------------

GROUP_INORGANIC = "inorganic"
GROUP_C1_C6 = "c1_c6"
GROUP_C7_PLUS = "c7_plus"

GROUP_TITLES = {
    GROUP_INORGANIC: "Inorganic (He, H2, N2, CO2, H2S)",
    GROUP_C1_C6: "Hydrocarbons C1-C6",
    GROUP_C7_PLUS: "C7+ fractions",
}


def component_catalog() -> list[dict]:
    """
    Все доступные компоненты БД в порядке `check_and_sort_composition`
    (сначала неорганика, затем углеводороды — оба по sequence_number).

    Returns
    -------
    list[dict]
        `{"name", "group", "sequence", "molar_mass", "gamma", "Tb"}` —
        M/gamma/Tb нужны как дефолты редактируемых колонок C7+ в форме.
    """
    db = _component_db()
    items: list[dict] = []
    for name in db["available_components"]:
        carbon = bool(db.get("carbon_flag", {}).get(name))
        c7 = bool(db.get("c7_plus_flag", {}).get(name))
        group = GROUP_INORGANIC if not carbon else (GROUP_C7_PLUS if c7 else GROUP_C1_C6)
        items.append({
            "name": name,
            "group": group,
            "sequence": db.get("sequence_number", {}).get(name, 0),
            "molar_mass": db.get("molar_mass", {}).get(name),
            "gamma": db.get("gamma", {}).get(name),
            "Tb": db.get("Tb", {}).get(name),
        })
    items.sort(key=lambda it: (0 if it["group"] == GROUP_INORGANIC else 1,
                               it["sequence"]))
    return items


# --- имена/валидация ---------------------------------------------------------

def _existing_ids(db_path: str) -> set:
    return set(read_model_store(db_path, missing_ok=True))


def suggest_model_id(name: str, db_path: str = "models.json") -> str:
    """
    Slug из имени модели: латиница/цифры/подчёркивания, верхний регистр
    (id должен быть валидным Python-идентификатором — это атрибут прокси
    `Composition.from_db`). Уникализируется суффиксом _2, _3, …
    """
    slug = re.sub(r"[^A-Za-z0-9]+", "_", (name or "").strip()).strip("_").upper()
    if not slug or not slug[0].isalpha():
        slug = "MODEL" + ("_" + slug if slug else "")
    existing = _existing_ids(db_path)
    candidate, n = slug, 2
    while candidate in existing:
        candidate = f"{slug}_{n}"
        n += 1
    return candidate


def _normalized_model_name(value: object) -> str:
    """Сравнивает отображаемые имена без различий регистра/лишних пробелов."""
    return " ".join(str(value or "").split()).casefold()


def _unique_copy_name(source_name: object, records: dict) -> str:
    """Возвращает удобное имя копии: ``Name (copy)``, затем ``(copy 2)``."""
    base = " ".join(str(source_name or "").split()) or "Unnamed model"
    existing = {_normalized_model_name(rec.get("Model_name"))
                for rec in records.values() if isinstance(rec, dict)}
    candidate = f"{base} (copy)"
    number = 2
    while _normalized_model_name(candidate) in existing:
        candidate = f"{base} (copy {number})"
        number += 1
    return candidate


def suggest_duplicate_name(model_id: str, db_path: str = "models.json") -> str:
    """Предлагает уникальное отображаемое имя копии модели."""
    records = read_model_store(db_path)
    source = records.get(model_id)
    if not isinstance(source, dict):
        raise KeyError(f"Model '{model_id}' not found in the database.")
    return _unique_copy_name(source.get("Model_name") or model_id, records)


def _unique_id_for_name(name: str, existing: set[str]) -> str:
    """Строит валидный уникальный model_id без повторного чтения файла."""
    slug = re.sub(r"[^A-Za-z0-9]+", "_", name.strip()).strip("_").upper()
    if not slug or not slug[0].isalpha():
        slug = "MODEL" + ("_" + slug if slug else "")
    candidate, number = slug, 2
    while candidate in existing:
        candidate = f"{slug}_{number}"
        number += 1
    return candidate


def duplicate_model(db_path: str, source_id: str, *, new_id: str | None = None,
                    new_name: str | None = None) -> str:
    """Создаёт независимую копию записи модели с новым id и именем.

    Копируются состав, EOS, пластовая температура, поле и корреляции. История
    сохранённых расчётов и GUI-workspace намеренно не копируются: после форка
    они могут относиться к исходной версии и вводить в заблуждение.
    Операция выполняется атомарно; перед записью повторно проверяются id и имя,
    поэтому параллельное создание не сможет затереть другую модель.
    """
    result_id: str | None = None

    def clone(records) -> None:
        nonlocal result_id
        source = records.get(source_id)
        if not isinstance(source, dict):
            raise KeyError(f"Model '{source_id}' not found in the database.")

        name = (new_name if new_name is not None
                else _unique_copy_name(source.get("Model_name") or source_id,
                                       records)).strip()
        if not name:
            raise ValueError("Model name is empty.")
        normalized = _normalized_model_name(name)
        for model_key, record in records.items():
            if (model_key != source_id and isinstance(record, dict)
                    and _normalized_model_name(record.get("Model_name")) == normalized):
                raise ValueError(f"Model name '{name}' already exists.")
            if model_key == source_id and _normalized_model_name(
                    record.get("Model_name")) == normalized:
                raise ValueError("Copy name must differ from the source model.")

        candidate_id = (new_id if new_id is not None
                        else _unique_id_for_name(name, set(records))).strip()
        if not candidate_id or not candidate_id.isidentifier():
            raise ValueError("Model id must be a valid identifier (letters/digits/_).")
        if candidate_id in records:
            raise ValueError(f"Model id '{candidate_id}' already exists.")

        now = datetime.now().isoformat(timespec="seconds")
        # Не тянем в память потенциально огромную историю результатов, которую
        # всё равно сбрасываем у новой модели.
        copy_record = deepcopy({key: value for key, value in source.items()
                                if key != "results"})
        copy_record["Model_name"] = name
        project_id = source.get("project_id")
        copy_record["project_id"] = (project_id.strip()
                                      if isinstance(project_id, str)
                                      and project_id.strip() else source_id)
        project_name = source.get("project_name")
        copy_record["project_name"] = (project_name if isinstance(project_name, str)
                                        and project_name.strip()
                                        else source.get("Model_name") or source_id)
        copy_record["created_at"] = now
        copy_record["updated_at"] = now
        copy_record["results"] = []
        records[candidate_id] = copy_record
        result_id = candidate_id

    update_model_store(db_path, clone)
    assert result_id is not None
    logger.info("Модель '%s' скопирована как '%s'", source_id, result_id)
    return result_id


def validate_new_model(model_id: str, name: str, zi: dict,
                       db_path: str = "models.json") -> list[str]:
    """Список ошибок формы (пустой список = всё ок). Тексты — на английском (UI)."""
    errors: list[str] = []
    if not (name or "").strip():
        errors.append("Model name is empty.")
    mid = (model_id or "").strip()
    if not mid or not mid.isidentifier():
        errors.append("Model id must be a valid identifier (letters/digits/_).")
    elif mid in _existing_ids(db_path):
        errors.append(f"Model id '{mid}' already exists in the database.")
    positive = {}
    for key, value in (zi or {}).items():
        try:
            number = float(value)
        except (TypeError, ValueError):
            errors.append(f"Mole fraction for '{key}' is not a number.")
            continue
        if not math.isfinite(number) or number < 0:
            errors.append(f"Mole fraction for '{key}' must be finite and >= 0.")
        elif number > 0:
            positive[key] = number
    if not positive:
        errors.append("Set at least one component fraction > 0.")
    available = set(_component_db()["available_components"])
    unknown = [k for k in positive if k not in available]
    if unknown:
        errors.append("Unknown components: " + ", ".join(unknown))
    return errors


# --- создание/сохранение -----------------------------------------------------

def create_model(db_path: str, model_id: str, name: str, field: Optional[str],
                 eos_value: str, t_res: float, zi: dict,
                 correlations: Optional[dict] = None,
                 c7_overrides: Optional[dict] = None) -> str:
    """
    Создаёт `Composition` из формы и сохраняет модель в models.json.

    Порядок важен: overrides M/gamma/Tb для C7+ применяются ДО
    `evaluate_composition_data` — эти свойства входят в корреляции C7+.

    Повреждённая или несовместимая база поднимает диагностируемую ошибку
    ModelStore и никогда не перезаписывается пустым содержимым.
    """
    zi_clean = {k: float(v) for k, v in zi.items() if v and float(v) > 0}
    composition = Composition(zi_clean, T_res=float(t_res),
                              eos_name=EOSType(eos_value))
    composition.normalize_composition()

    for comp_name, props in (c7_overrides or {}).items():
        if comp_name not in zi_clean:
            continue
        for prop, val in props.items():
            if val is not None:
                composition.edit_component_properties(comp_name, prop, float(val))

    selected_correlations = correlations or comp_svc.default_correlations()
    composition.evaluate_composition_data(selected_correlations)

    db = ModelJSONDB(db_path)
    db.export_and_save(model_id, name, composition.composition,
                       composition.composition_data, str(eos_value),
                       results=None, field=field or None, t_res=float(t_res),
                       correlations=selected_correlations,
                       project_id=model_id, project_name=name)
    logger.info("Модель '%s' создана и сохранена в %s", model_id, db_path)
    return model_id


# --- удаление модели ---------------------------------------------------------

def delete_model(db_path: str, model_id: str) -> bool:
    """
    Удаляет модель из models.json (с бэкапом `.bak`). Возвращает True, если
    модель была и удалена. Битый JSON пробрасывает исключение (не затираем).
    """
    if not Path(db_path).exists():
        return False

    class _ModelNotFound(Exception):
        pass

    def remove(data) -> None:
        if model_id not in data:
            raise _ModelNotFound
        del data[model_id]

    try:
        update_model_store(db_path, remove)
    except _ModelNotFound:
        return False
    logger.info("Модель '%s' удалена из %s", model_id, db_path)
    return True


def delete_project(db_path: str, project_id: str) -> bool:
    """Удаляет проект вместе со всеми его моделями.

    Для старых записей без ``project_id`` используется сам ключ модели как
    идентификатор проекта. Поэтому операция безопасна и для старого формата,
    где каждая модель фактически была отдельным проектом.
    """
    project_key = str(project_id).strip()
    if not project_key or not Path(db_path).exists():
        return False

    removed = 0

    class _ProjectNotFound(Exception):
        pass

    def remove(data) -> None:
        nonlocal removed
        for model_id, record in list(data.items()):
            record_project = (record.get("project_id") if isinstance(record, dict)
                              else None)
            effective_project = (record_project.strip()
                                 if isinstance(record_project, str)
                                 and record_project.strip() else model_id)
            if effective_project == project_key:
                del data[model_id]
                removed += 1
        if not removed:
            raise _ProjectNotFound

    try:
        update_model_store(db_path, remove)
    except _ProjectNotFound:
        return False
    if removed:
        logger.info("Проект '%s' удалён вместе с %d моделями из %s",
                    project_key, removed, db_path)
    return bool(removed)


# --- импорт Excel ------------------------------------------------------------

def excel_sheet_names(path: str) -> list[str]:
    """Имена листов Excel-файла (для выпадающего списка в предпросмотре)."""
    import pandas as pd
    return list(pd.ExcelFile(path).sheet_names)


def excel_preview(path: str, sheet: str, header: bool, max_rows: int = 25) -> dict:
    """
    Первые строки листа для предпросмотра (как их прочитает загрузчик).

    Returns
    -------
    dict
        `{"columns": [str, ...], "rows": [[str, ...], ...], "truncated": bool}`.
    """
    import pandas as pd
    hdr = 0 if header else None
    df = pd.read_excel(path, sheet_name=sheet, header=hdr, nrows=max_rows + 1)
    truncated = len(df) > max_rows
    df = df.head(max_rows)
    columns = [str(c) for c in df.columns]
    rows = [["" if pd.isna(v) else str(v) for v in rec]
            for rec in df.itertuples(index=False, name=None)]
    return {"columns": columns, "rows": rows, "truncated": truncated}


def import_excel(path: str, header: bool, sheet: str) -> dict:
    """
    Читает состав из Excel (`CompositionExcelLoader`) и сверяет имена с БД.

    Returns
    -------
    dict
        `{"recognized": {comp: zi}, "unrecognized": {name: value}}` —
        нераспознанные компоненты не теряются молча, их покажет форма.
    """
    raw = CompositionExcelLoader(path).load(header=header, sheet=sheet)
    available = set(_component_db()["available_components"])
    recognized = {k: float(v) for k, v in raw.items() if k in available}
    unrecognized = {k: v for k, v in raw.items() if k not in available}
    logger.info("Excel-импорт %s: %d распознано, %d нет",
                path, len(recognized), len(unrecognized))
    return {"recognized": recognized, "unrecognized": unrecognized}


# --- импорт E300 -------------------------------------------------------------

def _save_composition_as_model(db_path: str, model_id: str, name: str,
                               composition: Composition, eos_value: str,
                               t_res: float, field: Optional[str] = None,
                               correlations: Optional[dict] = None) -> None:
    """Сохраняет готовый ``Composition`` через единый безопасный ModelStore."""
    db = ModelJSONDB(db_path)
    db.export_and_save(model_id, name, composition.composition,
                       composition.composition_data, eos_value,
                       results=None, field=field, t_res=float(t_res),
                       correlations=correlations or comp_svc.default_correlations(),
                       project_id=model_id, project_name=name)


def preview_e300(path: str) -> dict:
    """Анализирует дек до сохранения и возвращает отчёт для подтверждения."""
    parsed = parse_e300(path)
    available = set(_component_db()["available_components"])
    recognized = [n for n in parsed["names"] if n in available]
    unrecognized = [n for n in parsed["names"] if n not in available]
    total = sum(float(v) for v in parsed["zi"].values())
    skipped = sum(float(parsed["zi"].get(n, 0.0)) for n in unrecognized)
    return {
        "name": Path(path).stem,
        "deck_eos": parsed["eos"],
        "recognized": recognized,
        "unrecognized": unrecognized,
        "total_zi": total,
        "skipped_zi": skipped,
        "warnings": list(parsed.get("warnings", [])),
    }


def import_e300(db_path: str, path: str, t_res: float = 373.15,
                *, allow_unrecognized: bool = False) -> dict:
    """
    Импортирует состав из E300-дека и сохраняет его отдельной моделью.

    Свойства из дека (MW/Tc/Pc/Vc/ω/shift/BIP) переносятся как есть: сначала
    `evaluate_composition_data` заполняет структуру, затем значения дека
    перезаписывают её (`edit_component_properties` при Tc/Pc/ω пересчитывает
    EOS-зависимые коэффициенты — модель остаётся согласованной). Модель
    создаётся на живом EOS Брусиловского независимо от ключа EOS в деке
    (в деке это лишь метка Eclipse — MPR/SRK).

    Returns
    -------
    dict
        `{"model_id", "name", "recognized": [...], "unrecognized": [...],
        "deck_eos": str|None}`.
    """
    parsed = parse_e300(path)
    available = set(_component_db()["available_components"])
    names = parsed["names"]
    recognized = [n for n in names if n in available]
    unrecognized = [n for n in names if n not in available]
    if unrecognized and not allow_unrecognized:
        raise ValueError(
            "E300 deck contains unrecognized components; preview and confirm "
            "before importing: " + ", ".join(unrecognized))

    zi_clean = {n: float(parsed["zi"][n]) for n in recognized
                if float(parsed["zi"].get(n, 0.0)) > 0}
    if not zi_clean:
        raise ValueError("E300 deck has no recognized components with a positive fraction")

    name = Path(path).stem
    model_id = suggest_model_id(name, db_path)

    composition = Composition(zi_clean, T_res=float(t_res), eos_name=EOSType.BRSEOS)
    composition.normalize_composition()
    composition.evaluate_composition_data(comp_svc.default_correlations())

    # перенос свойств из дека (только по распознанным компонентам)
    for data_key, col in parsed["props"].items():
        for comp_name, val in col.items():
            if comp_name in zi_clean and val is not None:
                try:
                    composition.edit_component_properties(comp_name, data_key, float(val))
                except KeyError:
                    logger.warning("E300-импорт: свойство %s недоступно", data_key)
    # перенос BIP
    bip = parsed["bip"]
    rec = list(zi_clean.keys())
    for ci in rec:
        for cj in rec:
            if ci != cj:
                composition.edit_bip_for_components(
                    ci, cj, float(bip.get(ci, {}).get(cj, 0.0)))

    _save_composition_as_model(db_path, model_id, name, composition,
                               str(EOSType.BRSEOS), t_res)
    logger.info("E300-импорт %s: модель '%s' (%d компонентов, %d нераспознано)",
                path, model_id, len(recognized), len(unrecognized))
    return {"model_id": model_id, "name": name, "recognized": recognized,
            "unrecognized": unrecognized, "deck_eos": parsed["eos"],
            "warnings": parsed.get("warnings", [])}


# --- сводка «что рассчитано» ------------------------------------------------

def calc_summary(summary, _legacy_session_ws: Optional[dict] = None) -> dict:
    """Счётчики сохранённых расчётов для колонки Calculated.

    Второй аргумент оставлен для совместимости: c R2.4 `gui_session.json`
    не участвует в представлении результатов.
    """
    records = [item for item in summary.results_brief if isinstance(item, dict)]
    computed = [item for item in records if item.get("has_result", True)]
    flashes = [item for item in computed
               if str(item.get("kind") or item.get("module") or "").lower()
               in {"flash", "flashresult"}]
    exps = [item for item in computed
            if str(item.get("kind") or "").lower() == "experiment"]
    envelopes = [item for item in computed
                 if str(item.get("kind") or "").lower() == "phase_envelope"]
    exp_kinds: dict[str, int] = {}
    for item in exps:
        kind = str(item.get("experiment_kind") or item.get("module")
                   or "experiment").upper()
        exp_kinds[kind] = exp_kinds.get(kind, 0) + 1
    return {
        "persisted": len(computed),
        "flashes": len(flashes),
        "experiments": len(exps),
        "envelopes": len(envelopes),
        "exp_kinds": exp_kinds,
        "stale": sum(1 for item in computed
                     if str(item.get("status") or "").upper() == "STALE"),
        "saved_at": getattr(summary, "workspace_saved_at", None),
    }
