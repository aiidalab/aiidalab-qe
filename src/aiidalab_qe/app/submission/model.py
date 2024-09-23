import traitlets as tl

from aiida import orm
from aiida.common import NotExistent
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.widgets import PwCodeResourceSetupWidget

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmissionModel(tl.HasTraits):
    input_structure = tl.Union(
        [
            tl.Instance(orm.StructureData),
            tl.Instance(HubbardStructureData),
        ],
        allow_none=True,
    )
    input_parameters = tl.Dict()

    process = tl.Instance(orm.WorkChainNode, allow_none=True)
    process_label = tl.Unicode("")
    process_description = tl.Unicode("")

    submission_blocker_messages = tl.Unicode("")

    installing_sssp = tl.Bool(allow_none=True)
    sssp_installed = tl.Bool(allow_none=True)
    setting_up_qe = tl.Bool(allow_none=True)

    codes = tl.Dict(
        key_trait=tl.Unicode(),
        value_trait=tl.Instance(
            PwCodeResourceSetupWidget,
            default_value={},
        ),
        default_value={},
    )

    def set_model_state(self, parameters):
        if "resources" in parameters:
            parameters["codes"] = {
                key: {"code": value} for key, value in parameters["codes"].items()
            }
            parameters["codes"]["pw"]["nodes"] = parameters["resources"]["num_machines"]
            parameters["codes"]["pw"]["cpus"] = parameters["resources"][
                "num_mpiprocs_per_machine"
            ]
            parameters["codes"]["pw"]["parallelization"] = {
                "npool": parameters["resources"]["npools"]
            }
        self.set_selected_codes(parameters["codes"])
        if self.process:
            self.process_label = self.process.label
            self.process_description = self.process.description

    def add_code(self, name, code):
        self.codes[name] = code  # type: ignore

    def get_code(self, name):
        return self.codes.get(name)

    def get_selected_codes(self):
        return {
            name: code.parameters
            for name, code in self.codes.items()
            if code.layout.display != "none"  # TODO do this differently
        }

    def set_selected_codes(self, code_data):
        def get_code_uuid(code):
            if code is not None:
                try:
                    return orm.load_code(code).uuid
                except NotExistent:
                    return None

        with self.hold_trait_notifications():
            for name, code in self.codes.items():
                if name not in code_data:
                    continue
                code_options = [
                    option[1]
                    for option in code.code_selection.code_select_dropdown.options
                ]
                parameters = code_data[name]
                code_uuid = get_code_uuid(parameters["code"])
                if code_uuid in code_options:
                    parameters["code"] = code_uuid
                    code.parameters = parameters

    def reset(self):
        with self.hold_trait_notifications():
            self.input_structure = None  # TODO why?
            self.process = None
            self.process_label = ""
            self.process_description = ""
            self.submission_blocker_messages = ""
