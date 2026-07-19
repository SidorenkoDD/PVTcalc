"""Общие fixtures для тестов PVTcalc."""

from pathlib import Path

import pytest

from calc_core.Composition.Composition import Composition

# Абсолютный путь от расположения самого conftest, а не относительный от cwd:
# иначе fixtures резолвятся только при запуске pytest из корня репозитория
# (`pytest tests/` из другой директории молча падал бы на поиске файла).
MODELS_JSON_PATH = str(Path(__file__).resolve().parents[1] / "models.json")


@pytest.fixture
def models_db():
    """Снэпшот моделей флюидов из models.json (кэшируется pytest на время теста)."""
    return Composition.from_db(MODELS_JSON_PATH)


@pytest.fixture
def krsnl_composition(models_db):
    """Состав KRSNL_PVTSIM — 40 компонентов, единственный, для которого название
    намекает на историческую сверку с PVTSim (само сравнение в репозитории не зафиксировано)."""
    comp = models_db.KRSNL_PVTSIM
    comp.normalize_composition()
    return comp


@pytest.fixture
def przlm_composition(models_db):
    """Состав PRRZLM_MDT_TEST — 35 компонентов, второй сохранённый снэпшот."""
    comp = models_db.PRRZLM_MDT_TEST
    comp.normalize_composition()
    return comp
