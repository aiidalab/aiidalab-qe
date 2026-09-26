import contextlib
import logging
import typing as t

import traitlets as tl

from aiidalab_qe.common.mixins import HasModels, HasProcess
from aiidalab_qe.common.mvc import Model


@contextlib.contextmanager
def capture_control_errors(
    logger: logging.Logger,
) -> t.Generator[list[str], None, None]:
    """Collect error messages logged by AiiDA process-control functions.

    AiiDA control functions can log failures such as unreachable processes and
    return normally, so checking only for raised exceptions can mistake failure
    for success.
    """
    messages = []

    class ErrorCaptureHandler(logging.Handler):
        def __init__(self):
            super().__init__(level=logging.ERROR)

        def emit(self, record: logging.LogRecord):
            messages.append(record.getMessage())

    handler = ErrorCaptureHandler()
    logger.addHandler(handler)
    try:
        yield messages
    finally:
        logger.removeHandler(handler)


class ResultsSubModel(Model, HasProcess):
    """Base class for result sub-models that have process-related models."""


PM = t.TypeVar("PM", bound=HasProcess)


class HasProcessModels(HasModels[PM]):
    """Mixin class for models that have process-related models."""

    def _link_model(self, model: PM):
        super()._link_model(model)
        tl.dlink(
            (self, "process_uuid"),
            (model, "process_uuid"),
        )
        tl.dlink(
            (self, "monitor_counter"),
            (model, "monitor_counter"),
        )
        tl.dlink(
            (self, "daemon_is_running"),
            (model, "daemon_is_running"),
        )
        tl.dlink(
            (self, "daemon_status_known"),
            (model, "daemon_status_known"),
        )
