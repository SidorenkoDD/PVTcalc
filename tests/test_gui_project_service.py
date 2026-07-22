"""
Тесты сервиса страницы Projects: каталог компонентов, создание модели
(roundtrip через Composition.from_db, включая новое поле T_res), защита от
затирания базы, импорт Excel. Без DearPyGui — идут в CI.
"""

import json
from pathlib import Path

import pandas as pd
import pytest

from calc_core.Composition.Composition import Composition
from calc_core.Utils.ModelStore import ModelStoreCorruptError
from gui.app_state import AppState
from gui.services import project_service as svc
from gui.services.model_repository import ModelRepository
from gui.workspace_codec import restore_workspace

MODELS_JSON = Path(__file__).resolve().parent / "fixtures" / "models.json"

_ZI = {"N2": 0.01, "CO2": 0.02, "C1": 0.60, "C2": 0.08, "C3": 0.05,
       "C6": 0.04, "C7": 0.20}


# --- каталог компонентов -----------------------------------------------------

def test_component_catalog_full_and_grouped():
    cat = svc.component_catalog()
    names = [c["name"] for c in cat]
    assert len(names) == 52  # 5 неорганических + C1..C6 (8) + C7..C45 (39)
    assert names[0] == "He"                       # неорганика первой
    groups = {c["name"]: c["group"] for c in cat}
    assert groups["N2"] == svc.GROUP_INORGANIC
    assert groups["C3"] == svc.GROUP_C1_C6
    assert groups["C7"] == svc.GROUP_C7_PLUS
    c7 = next(c for c in cat if c["name"] == "C7")
    assert c7["molar_mass"] and c7["gamma"] and c7["Tb"]  # дефолты для формы


# --- имена/валидация ---------------------------------------------------------

def test_suggest_model_id_slug_and_unique(tmp_path):
    db = tmp_path / "m.json"
    db.write_text(
        '{"MY_FLUID": {"composition": {"C1": 1.0}, '
        '"composition_data": {}}}', encoding="utf-8")
    assert svc.suggest_model_id("My fluid!", str(db)) == "MY_FLUID_2"
    assert svc.suggest_model_id("2 phases", str(db)).startswith("MODEL")
    assert svc.suggest_model_id("нефть", str(db)).startswith("MODEL")


def test_duplicate_model_uses_unique_name_and_resets_results(tmp_path):
    db = str(tmp_path / "models.json")
    svc.create_model(db, "A", "Original", "Field", "PREOS", 373.15, dict(_ZI))
    svc.create_model(db, "B", "Original (copy)", None, "PREOS", 373.15, dict(_ZI))

    from calc_core.Utils.ModelStore import update_model_store

    def add_result(records):
        records["A"]["results"] = [{"module": "Flash", "data": {"x": 1}}]

    update_model_store(db, add_result)
    copy_id = svc.duplicate_model(db, "A")
    data = json.loads(Path(db).read_text(encoding="utf-8"))
    assert copy_id in data and copy_id not in {"A", "B"}
    assert data[copy_id]["Model_name"] == "Original (copy 2)"
    assert data[copy_id]["composition"] == data["A"]["composition"]
    assert data[copy_id]["correlations"] == data["A"]["correlations"]
    assert data[copy_id]["project_id"] == "A"
    assert data[copy_id]["project_name"] == "Original"
    assert data[copy_id]["results"] == []
    assert data[copy_id]["created_at"] != data["A"]["created_at"]


def test_duplicate_model_rejects_existing_id_or_name(tmp_path):
    db = str(tmp_path / "models.json")
    svc.create_model(db, "A", "Original", None, "PREOS", 373.15, dict(_ZI))
    with pytest.raises(ValueError, match="name"):
        svc.duplicate_model(db, "A", new_id="C", new_name=" original ")
    with pytest.raises(ValueError, match="id"):
        svc.duplicate_model(db, "A", new_id="A", new_name="Another")


def test_validate_new_model_errors(tmp_path):
    db = tmp_path / "m.json"
    db.write_text(
        '{"X": {"composition": {"C1": 1.0}, '
        '"composition_data": {}}}', encoding="utf-8")
    errors = svc.validate_new_model("X", "", {}, str(db))
    joined = " ".join(errors)
    assert "already exists" in joined
    assert "name is empty" in joined
    assert "at least one component" in joined
    assert svc.validate_new_model("OK_ID", "Name", {"C1": 1.0}, str(db)) == []
    assert "Unknown components" in " ".join(
        svc.validate_new_model("OK2", "Name", {"C1+": 1.0}, str(db)))


# --- создание модели ---------------------------------------------------------

def test_create_model_roundtrip(tmp_path):
    db = str(tmp_path / "models.json")
    svc.create_model(db, "NEW_FLUID", "New fluid", "TestField", "PREOS",
                     t_res=380.0, zi=dict(_ZI))
    # roundtrip через штатный from_db
    proxy = Composition.from_db(db)
    assert proxy.list_models() == ["NEW_FLUID"]
    comp = proxy.NEW_FLUID
    assert comp.T == pytest.approx(380.0)          # T_res сохранился (новое поле)
    assert comp.eos_name.value == "PREOS"
    assert comp.composition["C1"] == pytest.approx(0.60, abs=1e-6)
    # свойства посчитаны (C7 через корреляции)
    assert comp.composition_data["critical_temperature"]["C7"] > 0
    # сводка репозитория видит новое поле
    s = ModelRepository(db_path=db).list_models()[0]
    assert s.t_res == pytest.approx(380.0)
    assert s.title == "New fluid" and s.field_name == "TestField"


def test_create_model_drops_zero_zi_and_applies_overrides(tmp_path):
    db = str(tmp_path / "models.json")
    zi = dict(_ZI)
    zi["H2S"] = 0.0                                 # нулевая доля отбрасывается
    svc.create_model(db, "F2", "F2", None, "PREOS", 373.15, zi,
                     c7_overrides={"C7": {"molar_mass": 99.0}})
    comp = Composition.from_db(db).F2
    assert "H2S" not in comp.composition
    assert comp.composition_data["molar_mass"]["C7"] == pytest.approx(99.0)


def test_create_model_appends_not_overwrites(tmp_path):
    db = str(tmp_path / "models.json")
    svc.create_model(db, "A", "A", None, "PREOS", 373.15, dict(_ZI))
    svc.create_model(db, "B", "B", None, "PREOS", 373.15, dict(_ZI))
    data = json.loads(Path(db).read_text(encoding="utf-8"))
    assert set(data.keys()) == {"A", "B"}
    assert Path(db + ".bak").exists()               # бэкап перед второй записью


def test_create_model_refuses_on_corrupt_db(tmp_path):
    db = tmp_path / "models.json"
    db.write_text("{corrupted", encoding="utf-8")   # непустой, но битый JSON
    with pytest.raises(ModelStoreCorruptError, match="строка 1, столбец 2"):
        svc.create_model(str(db), "X", "X", None, "PREOS", 373.15, dict(_ZI))
    assert db.read_text(encoding="utf-8") == "{corrupted"  # база не тронута


# --- импорт Excel ------------------------------------------------------------

def _write_xlsx(path, rows, header):
    df = pd.DataFrame(rows, columns=["component", "zi"] if header else None)
    df.to_excel(path, index=False, header=header)


def test_import_excel_with_header(tmp_path):
    x = tmp_path / "comp.xlsx"
    _write_xlsx(x, [["C1", 0.7], ["C7", 0.2], ["Хитрый+", 0.1]], header=True)
    res = svc.import_excel(str(x), header=True, sheet="Sheet1")
    assert res["recognized"] == {"C1": 0.7, "C7": 0.2}
    assert "Хитрый+" in res["unrecognized"]         # не потеряли молча


def test_import_excel_no_header(tmp_path):
    x = tmp_path / "comp.xlsx"
    _write_xlsx(x, [["C1", 0.9], ["CO2", 0.1]], header=False)
    res = svc.import_excel(str(x), header=False, sheet="Sheet1")
    assert res["recognized"] == {"C1": 0.9, "CO2": 0.1}
    assert res["unrecognized"] == {}


def test_import_excel_no_header_honors_selected_sheet(tmp_path):
    x = tmp_path / "multi-no-header.xlsx"
    with pd.ExcelWriter(x) as writer:
        pd.DataFrame([["C1", 1.0]]).to_excel(
            writer, sheet_name="first", index=False, header=False)
        pd.DataFrame([["CO2", 1.0]]).to_excel(
            writer, sheet_name="selected", index=False, header=False)
    res = svc.import_excel(str(x), header=False, sheet="selected")
    assert res["recognized"] == {"CO2": 1.0}


def test_import_excel_bad_values_raise(tmp_path):
    x = tmp_path / "bad.xlsx"
    _write_xlsx(x, [["C1", "not-a-number"]], header=True)
    with pytest.raises(Exception):
        svc.import_excel(str(x), header=True, sheet="Sheet1")


def _write_multi_sheet(path):
    import openpyxl
    wb = openpyxl.Workbook()
    ws1 = wb.active
    ws1.title = "compo"
    ws1.append(["component", "zi"])
    ws1.append(["C1", 0.7])
    ws1.append(["C7", 0.3])
    ws2 = wb.create_sheet("other")
    ws2.append(["a", "b"])
    ws2.append(["x", 1])
    wb.save(path)


def test_excel_sheet_names_and_preview(tmp_path):
    x = tmp_path / "multi.xlsx"
    _write_multi_sheet(x)
    assert svc.excel_sheet_names(str(x)) == ["compo", "other"]

    prev = svc.excel_preview(str(x), "compo", header=True)
    assert prev["columns"] == ["component", "zi"]
    assert prev["rows"][0] == ["C1", "0.7"]
    assert prev["truncated"] is False
    # выбор листа реально применяется
    assert svc.excel_preview(str(x), "other", header=True)["columns"] == ["a", "b"]


def test_excel_preview_truncates(tmp_path):
    x = tmp_path / "big.xlsx"
    pd.DataFrame([["C1", i] for i in range(40)],
                 columns=["c", "v"]).to_excel(x, index=False)
    prev = svc.excel_preview(str(x), "Sheet1", header=True, max_rows=10)
    assert len(prev["rows"]) == 10
    assert prev["truncated"] is True


# --- удаление модели ---------------------------------------------------------

def test_delete_model(tmp_path):
    db = str(tmp_path / "models.json")
    svc.create_model(db, "A", "A", None, "PREOS", 373.15, dict(_ZI))
    svc.create_model(db, "B", "B", None, "PREOS", 373.15, dict(_ZI))
    assert svc.delete_model(db, "A") is True
    assert list(json.loads(Path(db).read_text(encoding="utf-8"))) == ["B"]
    assert Path(db + ".bak").exists()
    # повторное удаление / отсутствующая модель / нет файла
    assert svc.delete_model(db, "A") is False
    assert svc.delete_model(str(tmp_path / "нет.json"), "A") is False


def test_delete_project_removes_all_model_variants(tmp_path):
    db = str(tmp_path / "models.json")
    svc.create_model(db, "A", "Original", None, "PREOS", 373.15, dict(_ZI))
    copy_id = svc.duplicate_model(db, "A", new_id="A_COPY",
                                  new_name="Original (copy)")
    svc.create_model(db, "OTHER", "Other", None, "PREOS", 373.15, dict(_ZI))

    assert svc.delete_project(db, "A") is True
    data = json.loads(Path(db).read_text(encoding="utf-8"))
    assert set(data) == {"OTHER"}
    assert svc.delete_project(db, "A") is False
    assert copy_id == "A_COPY"


def test_repository_saves_edited_composition_and_preserves_metadata(tmp_path):
    db = str(tmp_path / "models.json")
    svc.create_model(db, "A", "Original", "Field", "PREOS", 373.15, dict(_ZI))
    repo = ModelRepository(db)
    comp = repo.load_composition("A")
    comp.composition["C1"] = 0.61
    comp.normalize_composition()
    repo.save_composition("A", comp, {"critical_temperature": "pedersen"})
    raw = json.loads(Path(db).read_text(encoding="utf-8"))["A"]
    assert raw["Model_name"] == "Original" and raw["Field"] == "Field"
    assert raw["updated_at"]
    assert raw["correlations"]["critical_temperature"] == "pedersen"
    assert Path(db + ".bak").exists()


def test_save_model_persists_workspace_results_and_provenance(tmp_path):
    db = str(tmp_path / "models.json")
    svc.create_model(db, "A", "Original", "Field", "PREOS", 373.15, dict(_ZI))
    state = AppState(ModelRepository(db))
    state.refresh_model_list()
    state.enter_model("A")
    node_id = state.new_envelope({"method": "ssm", "t_min_c": 10.0})
    assert node_id is not None
    state.set_node_result(node_id, {"points": [[10.0, 42.0]]})
    state.save_model()

    raw = json.loads(Path(db).read_text(encoding="utf-8"))["A"]
    stored = raw["workspace"]
    assert stored["schema_version"] == 1
    assert stored["provenance"]["composition_sha256"]
    assert stored["snapshot"]["nodes"][0]["node_id"] == node_id
    assert raw["results"][0]["has_result"] is True

    fresh = AppState(ModelRepository(db))
    fresh.refresh_model_list()
    model = fresh.ensure_model_loaded("A")
    saved_workspace = fresh.load_saved_workspace("A")
    assert saved_workspace is not None
    restore_workspace(model.variants["base"], saved_workspace)
    restored = model.variants["base"].nodes[node_id]
    assert restored.result == {"points": [[10.0, 42.0]]}
    summary = fresh.models["A"].summary
    assert summary is not None
    assert summary.workspace_saved_at == stored["saved_at"]
    info = svc.calc_summary(summary, {"flashes": [{"result": {"temporary": True}}]})
    assert info["persisted"] == 1 and info["flashes"] == 0
    assert info["envelopes"] == 1
