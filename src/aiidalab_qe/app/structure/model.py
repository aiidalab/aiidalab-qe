import traitlets as tl

from aiida import orm
from aiidalab_qe.common.mixins import Confirmable, HasTraitsAndMixins


class StructureModel(HasTraitsAndMixins, Confirmable):
    structure = tl.Instance(
        orm.StructureData,
        allow_none=True,
    )
    structure_name = tl.Unicode("")
    manager_output = tl.Unicode("")
    message_area = tl.Unicode("")

    def update_widget_text(self):
        if self.structure is None:
            self.structure_name = ""
            self.message_area = ""
        else:
            self.manager_output = ""
            self.structure_name = str(self.structure.get_formula())

    def reset(self):
        self.structure = None
        self.structure_name = ""
        self.manager_output = ""
        self.message_area = ""
