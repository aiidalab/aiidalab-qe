import traitlets as tl

from aiida import orm
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData


class StructureModel(tl.HasTraits):
    structure = tl.Instance(
        orm.StructureData,
        allow_none=True,
    )
    confirmed_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )

    confirmed = tl.Bool(False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.observe(
            self._unconfirm,
            "confirmed_structure",
        )

    def reset(self):
        self.confirmed_structure = None

    def _unconfirm(self, _):
        self.confirmed = False


struct_model = StructureModel()
