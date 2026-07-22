"""Фоновые задачи адресуются полным NodeRef и не зависят от DPG."""

import threading
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


def test_cancelled_operation_receives_token_and_progress():
    jobs = CalculationCoordinator()
    ref = NodeRef("model", "base", "ssm_1")
    started = threading.Event()

    def operation(token, progress):
        started.set()
        progress(0.25, "working")
        while True:
            token.throw_if_cancelled()
            time.sleep(0.001)

    job = jobs.start(ref, "envelope", operation)
    assert started.wait(1.0)
    assert job.progress == 0.25
    assert job.progress_message == "working"
    assert jobs.cancel() is True
    assert job.done.wait(1.0)

    finished = jobs.take_finished()
    assert finished is not None
    assert finished.cancelled is True
    assert finished.result is None
    assert finished.error is None
