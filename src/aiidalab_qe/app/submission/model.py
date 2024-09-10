import traitlets

from aiida import orm
from aiidalab_widgets_base import WizardAppWidgetStep


class SubmissionModel(traitlets.HasTraits):
    state = traitlets.UseEnum(WizardAppWidgetStep.State)
    prev_step_state = traitlets.UseEnum(WizardAppWidgetStep.State)

    input_structure = traitlets.Instance(orm.StructureData, allow_none=True)
    input_parameters = traitlets.Dict()
    process = traitlets.Instance(orm.WorkChainNode, allow_none=True)


submit_model = SubmissionModel()
