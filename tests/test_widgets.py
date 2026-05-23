import asyncio
from threading import Thread

from aiidalab_qe.common import widgets
from aiidalab_qe.common.widgets import CalcJobOutputFollower


class _SealedCalcJob:
    is_sealed = True


def test_calcjob_output_follower_sets_event_loop_in_worker_thread(monkeypatch):
    follower = CalcJobOutputFollower()
    errors = []

    monkeypatch.setattr(widgets, "load_node", lambda _uuid: _SealedCalcJob())

    def fetch_output(_calcjob):
        try:
            asyncio.get_event_loop()
        except RuntimeError as error:
            errors.append(error)
            raise
        return ["output line"]

    monkeypatch.setattr(follower, "_fetch_output", fetch_output)

    thread = Thread(target=follower._push_output, args=("uuid", 0))
    thread.start()
    thread.join(timeout=5)

    assert not thread.is_alive()
    assert errors == []
    assert follower._output_queue.get_nowait() == ["output line"]
    assert follower._output_queue.get_nowait() is follower._EOF
