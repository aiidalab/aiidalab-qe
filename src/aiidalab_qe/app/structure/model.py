import traitlets as tl

from aiida import orm
from aiidalab_widgets_base import WizardAppWidgetStep


class StructureModel(tl.HasTraits):
    state = tl.UseEnum(
        enum_class=WizardAppWidgetStep.State,
        default_value=WizardAppWidgetStep.State.INIT,
    )

    confirmed_structure = tl.Instance(orm.StructureData, allow_none=True)


struct_model = StructureModel()
