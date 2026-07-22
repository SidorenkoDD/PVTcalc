import pytest

from gui.services import lab_data_service as lab_svc


def test_project_catalog_shares_data_and_keeps_model_scope_private(tmp_path):
    db_path = tmp_path / "models.json"
    columns = ["pressure", "Bo"]
    shared = lab_svc.create_dataset(
        db_path, "project_a", title="DLE reference", experiment_kind="dle",
        columns=columns, rows=[[300.0, 1.2]], scope="project",
    )
    local = lab_svc.create_dataset(
        db_path, "project_a", title="Tuned DLE", experiment_kind="dle",
        columns=columns, rows=[[300.0, 1.3]], scope="model", model_id="tuned",
    )

    assert lab_svc.catalog_path(db_path).exists()
    assert [item["dataset_id"] for item in lab_svc.list_datasets(
        db_path, "project_a", experiment_kind="dle")] == [shared["dataset_id"]]
    assert [item["dataset_id"] for item in lab_svc.list_datasets(
        db_path, "project_a", experiment_kind="dle", model_id="tuned")] == [
            shared["dataset_id"], local["dataset_id"],
        ]
    assert lab_svc.effective_lab_data(
        db_path, "project_a", "other", {"lab_data_ref": local["dataset_id"],
                                            "lab_data": {"rows": [[1.0]]}},
    ) == {"rows": [[1.0]]}
    resolved = lab_svc.effective_lab_data(
        db_path, "project_a", "tuned", {"lab_data_ref": local["dataset_id"]},
    )
    assert resolved is not None
    assert resolved["rows"] == [[300.0, 1.3]]


def test_model_and_project_catalog_cleanup(tmp_path):
    db_path = tmp_path / "models.json"
    common = dict(experiment_kind="dle", columns=["pressure"], rows=[[200.0]])
    lab_svc.create_dataset(db_path, "project_a", title="shared", **common)
    lab_svc.create_dataset(db_path, "project_a", title="model", scope="model",
                           model_id="model_a", **common)

    lab_svc.delete_model_datasets(db_path, "project_a", "model_a")
    assert len(lab_svc.list_datasets(db_path, "project_a")) == 1
    lab_svc.delete_project_datasets(db_path, "project_a")
    assert lab_svc.list_datasets(db_path, "project_a") == []


def test_dataset_update_preserves_scope_and_manual_conditions(tmp_path):
    db_path = tmp_path / "models.json"
    dataset = lab_svc.create_dataset(
        db_path, "project_a", title="DLE initial", experiment_kind="dle",
        columns=["pressure", "Bo"], rows=[[300.0, 1.2]],
        conditions={"T_c": 80.0, "P_res": 350.0},
    )

    updated = lab_svc.update_dataset(
        db_path, "project_a", dataset["dataset_id"], title="DLE revised",
        columns=["pressure", "Bo"], rows=[[250.0, 1.1]],
        conditions={"T_c": 85.0, "P_res": 300.0},
    )

    assert updated is not None
    assert updated["title"] == "DLE revised"
    assert updated["conditions"] == {"T_c": 85.0, "P_res": 300.0}
    assert updated["rows"] == [[250.0, 1.1]]


def test_dataset_can_be_deleted_within_its_scope(tmp_path):
    db_path = tmp_path / "models.json"
    dataset = lab_svc.create_dataset(
        db_path, "project_a", title="Private DLE", experiment_kind="dle",
        columns=["pressure"], rows=[[200.0]], scope="model", model_id="model_a",
    )

    assert lab_svc.delete_dataset(
        db_path, "project_a", dataset["dataset_id"], model_id="model_a") is True
    assert lab_svc.list_datasets(db_path, "project_a", model_id="model_a") == []


def test_corrupt_catalog_is_not_overwritten(tmp_path):
    db_path = tmp_path / "models.json"
    path = lab_svc.catalog_path(db_path)
    path.write_text("{broken", encoding="utf-8")

    with pytest.raises(lab_svc.LabDataStoreError):
        lab_svc.create_dataset(
            db_path, "project_a", title="DLE", experiment_kind="dle",
            columns=["pressure"], rows=[[100.0]],
        )

    assert path.read_text(encoding="utf-8") == "{broken"
