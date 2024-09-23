import traitlets as tl

from aiida import orm
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.common.widgets import PwCodeResourceSetupWidget


class SubmissionModel(tl.HasTraits):
    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )

    process = tl.Instance(orm.WorkChainNode, allow_none=True)
    input_parameters = tl.Dict()

    codes = tl.Dict(
        key_trait=tl.Unicode(),
        value_trait=tl.Dict(
            key_trait=tl.Unicode(),
            value_trait=tl.Instance(
                PwCodeResourceSetupWidget,
                default_value={},
            ),
            default_value={},
        ),
        default_value={},
    )

    def get_model_state(self):
        return {}

    def set_model_state(self, parameters):
        pass

    def reset(self):
        self.input_parameters = {}
        self.codes = {}
