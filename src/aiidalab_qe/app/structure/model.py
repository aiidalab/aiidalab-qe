import traitlets as tl

from aiidalab_qe.common.mixins import Confirmable, HasInputStructure, HasTraitsAndMixins


class StructureModel(
    HasTraitsAndMixins,
    HasInputStructure,
    Confirmable,
):
    structure_name = tl.Unicode("")
    manager_output = tl.Unicode("")
    message_area = tl.Unicode("")

    def update_widget_text(self):
        if not self.has_structure:
            self.structure_name = ""
            self.message_area = ""
        else:
            self.manager_output = ""
            self.structure_name = str(self.input_structure.get_formula())

    def reset(self):
        self.input_structure = None
        self.structure_name = ""
        self.manager_output = ""
        self.message_area = ""
