import typing as t

import traitlets as tl

from aiidalab_qe.common.mixins import HasModels, HasProcess
from aiidalab_qe.common.mvc import Model


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
