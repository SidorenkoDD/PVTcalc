"""Фоновые задачи адресуются полным NodeRef и не зависят от DPG."""

import time

import pytest

from gui.app_state import NodeRef
from gui.calculation_coordinator import CalculationCoordinator


def test_job_completes_and_is_taken():
    jobs = CalculationCoordinator()
    ref = NodeRef("model", "base", "flash_1")
    job = jobs.start(ref, "flash", lambda: 42)
    assert jobs.busy
    assert job.done.wait(1.0)
    finished = jobs.take_finished()
    assert finished.node_ref == ref
    assert finished.result == 42
    assert not jobs.busy


def test_only_one_job_and_cancel_discards_flag():
    jobs = CalculationCoordinator()
    ref = NodeRef("model", "base", "env_1")
    job = jobs.start(ref, "envelope", lambda: (time.sleep(0.02), "done")[1])
    with pytest.raises(RuntimeError):
        jobs.start(ref, "second", lambda: None)
    assert jobs.cancel() is True
    assert job.done.wait(1.0)
    assert jobs.take_finished().cancelled is True
