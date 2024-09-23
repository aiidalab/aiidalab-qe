import traitlets as tl

from aiida import orm
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData


class StructureModel(tl.HasTraits):
    confirmed_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )

    def reset(self):
        self.confirmed_structure = None

    @property
    def is_confirmed(self):
        return self.confirmed_structure is not None


struct_model = StructureModel()
