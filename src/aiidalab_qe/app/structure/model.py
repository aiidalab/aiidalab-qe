import traitlets as tl

from aiida import orm
from aiidalab_widgets_base import WizardAppWidgetStep


class StructureModel(tl.HasTraits):
    state = tl.UseEnum(WizardAppWidgetStep.State)

    confirmed_structure = tl.Instance(orm.StructureData, allow_none=True)


struct_model = StructureModel()
