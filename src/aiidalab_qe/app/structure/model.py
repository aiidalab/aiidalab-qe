import traitlets as tl

from aiida import orm


class StructureModel(tl.HasTraits):
    structure = tl.Instance(
        orm.StructureData,
        allow_none=True,
    )
    structure_name = tl.Unicode("")
    manager_output = tl.Unicode("")
    message_area = tl.Unicode("")
    confirmed = tl.Bool(False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.observe(
            self._unconfirm,
            "structure",
        )

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

    def _unconfirm(self, _):
        self.confirmed = False
