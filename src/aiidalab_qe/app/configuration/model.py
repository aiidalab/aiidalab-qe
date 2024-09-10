import traitlets

from aiida import orm
from aiidalab_widgets_base import WizardAppWidgetStep


class ConfigurationModel(traitlets.HasTraits):
    state = traitlets.UseEnum(WizardAppWidgetStep.State)
    prev_step_state = traitlets.UseEnum(WizardAppWidgetStep.State)

    input_structure = traitlets.Instance(orm.StructureData, allow_none=True)
    configuration_parameters = traitlets.Dict()
    workchain_protocol = traitlets.Unicode()
    spin_type = traitlets.Unicode()
    electronic_type = traitlets.Unicode()


config_model = ConfigurationModel()
