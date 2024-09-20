import traitlets as tl

from aiida import orm
from aiidalab_widgets_base import WizardAppWidgetStep


class SubmissionModel(tl.HasTraits):
    state = tl.UseEnum(
        enum_class=WizardAppWidgetStep.State,
        default_value=WizardAppWidgetStep.State.INIT,
    )
    previous_step_state = tl.UseEnum(WizardAppWidgetStep.State)

    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    input_parameters = tl.Dict()
    process = tl.Instance(orm.WorkChainNode, allow_none=True)


submit_model = SubmissionModel()
