import traitlets

from aiida import orm
from aiidalab_widgets_base import WizardAppWidgetStep


class StructureModel(traitlets.HasTraits):
    state = traitlets.UseEnum(WizardAppWidgetStep.State)

    confirmed_structure = traitlets.Instance(orm.StructureData, allow_none=True)


struct_model = StructureModel()
