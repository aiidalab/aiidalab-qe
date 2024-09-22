import traitlets as tl

from aiida import orm


class StructureModel(tl.HasTraits):
    confirmed_structure = tl.Instance(orm.StructureData, allow_none=True)

    def reset(self):
        self.confirmed_structure = None

    @property
    def is_confirmed(self):
        return self.confirmed_structure is not None


struct_model = StructureModel()
