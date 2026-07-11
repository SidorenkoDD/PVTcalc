"""Общие fixtures для тестов PVTcalc."""

import pytest

from _src.Composition.Composition import Composition

MODELS_JSON_PATH = "models.json"


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
