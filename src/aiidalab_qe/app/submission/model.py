import traitlets as tl

from aiida import orm


class SubmissionModel(tl.HasTraits):
    input_structure = tl.Instance(orm.StructureData, allow_none=True)
    input_parameters = tl.Dict()
    process = tl.Instance(orm.WorkChainNode, allow_none=True)
    codes = tl.Dict()


submit_model = SubmissionModel()
