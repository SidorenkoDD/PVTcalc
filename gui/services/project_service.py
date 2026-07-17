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

import json
import logging
import re
import shutil
from pathlib import Path
from typing import Optional

from calc_core.Composition.Composition import Composition
from calc_core.EOS.BaseEOS import EOSType
from calc_core.Utils.CompositionLoader import CompositionExcelLoader
from calc_core.Utils.E300Import import parse_e300
from calc_core.Utils.Export import ModelJSONDB
from calc_core.Utils.JsonDBReader import JsonDBReader

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
    p = Path(db_path)
    if not p.exists():
        return set()
    try:
        with open(p, "r", encoding="utf-8") as f:
            return set(json.load(f).keys())
    except json.JSONDecodeError:
        return set()


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
    positive = {k: v for k, v in (zi or {}).items() if v and float(v) > 0}
    if not positive:
        errors.append("Set at least one component fraction > 0.")
    available = set(_component_db()["available_components"])
    unknown = [k for k in positive if k not in available]
    if unknown:
        errors.append("Unknown components: " + ", ".join(unknown))
    return errors


# --- создание/сохранение -----------------------------------------------------

def create_model(db_path: str, model_id: str, name: str, field: str,
                 eos_value: str, t_res: float, zi: dict,
                 correlations: Optional[dict] = None,
                 c7_overrides: Optional[dict] = None) -> str:
    """
    Создаёт `Composition` из формы и сохраняет модель в models.json.

    Порядок важен: overrides M/gamma/Tb для C7+ применяются ДО
    `evaluate_composition_data` — эти свойства входят в корреляции C7+.

    Raises
    ------
    RuntimeError
        Если существующий непустой models.json не удалось загрузить —
        прерываемся, чтобы `save()` не затёр базу.
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

    composition.evaluate_composition_data(
        correlations or comp_svc.default_correlations())

    # защита от затирания базы: файл есть и непуст, а load() дал пусто
    p = Path(db_path)
    db = ModelJSONDB(db_path)
    if p.exists() and p.stat().st_size > 2 and not db._db:
        raise RuntimeError(
            "models database exists but could not be loaded - aborting save "
            "to avoid data loss")
    if p.exists():
        shutil.copyfile(p, str(p) + ".bak")

    db.export_and_save(model_id, name, composition.composition,
                       composition.composition_data, str(eos_value),
                       results=None, field=field or None, t_res=float(t_res))
    logger.info("Модель '%s' создана и сохранена в %s", model_id, db_path)
    return model_id


# --- удаление модели ---------------------------------------------------------

def delete_model(db_path: str, model_id: str) -> bool:
    """
    Удаляет модель из models.json (с бэкапом `.bak`). Возвращает True, если
    модель была и удалена. Битый JSON пробрасывает исключение (не затираем).
    """
    p = Path(db_path)
    if not p.exists():
        return False
    with open(p, "r", encoding="utf-8") as f:
        data = json.load(f)
    if model_id not in data:
        return False
    shutil.copyfile(p, str(p) + ".bak")
    del data[model_id]
    with open(p, "w", encoding="utf-8") as f:
        json.dump(data, f, indent=2, ensure_ascii=False)
    logger.info("Модель '%s' удалена из %s", model_id, db_path)
    return True


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
                               t_res: float, field: Optional[str] = None) -> None:
    """Сохраняет готовый `Composition` как модель (guard базы + бэкап `.bak`)."""
    p = Path(db_path)
    db = ModelJSONDB(db_path)
    if p.exists() and p.stat().st_size > 2 and not db._db:
        raise RuntimeError(
            "models database exists but could not be loaded - aborting save "
            "to avoid data loss")
    if p.exists():
        shutil.copyfile(p, str(p) + ".bak")
    db.export_and_save(model_id, name, composition.composition,
                       composition.composition_data, eos_value,
                       results=None, field=field, t_res=float(t_res))


def import_e300(db_path: str, path: str, t_res: float = 373.15) -> dict:
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
            "unrecognized": unrecognized, "deck_eos": parsed["eos"]}


# --- сводка «что рассчитано» ------------------------------------------------

def calc_summary(summary, session_ws: Optional[dict]) -> dict:
    """
    Счётчики для колонки Calculated страницы Projects.

    Parameters
    ----------
    summary : ModelSummary
        Сводка модели (persisted-результаты из models.json).
    session_ws : dict | None
        Workspace модели из сессии v2 (`{flashes, experiments, ...}`).
    """
    ws = session_ws or {}
    flashes = [f for f in ws.get("flashes", []) if f.get("result")]
    exps = [e for e in ws.get("experiments", []) if e.get("result")]
    exp_kinds: dict[str, int] = {}
    for e in exps:
        k = str(e.get("kind", "?")).upper()
        exp_kinds[k] = exp_kinds.get(k, 0) + 1
    return {
        "persisted": len(summary.results_brief),
        "flashes": len(flashes),
        "experiments": len(exps),
        "exp_kinds": exp_kinds,
    }
