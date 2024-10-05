from __future__ import annotations

import typing as t

import traitlets as tl

from aiida import orm
from aiida.common import NotExistent
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class CodeModel(tl.HasTraits):
    is_active = tl.Bool(False)
    selected = tl.Bool(False)
    options = tl.List(
        trait=tl.Tuple(tl.Unicode, tl.Unicode),  # code option (label, uuid)
        default_value=[],
    )
    parameters = tl.Dict()

    def __init__(
        self,
        name,
        description,
        default_calc_job_plugin,
        setup_widget_class=QEAppComputationalResourcesWidget,
    ):
        self.name = name
        self.description = description
        self.default_calc_job_plugin = default_calc_job_plugin
        self.setup_widget_class = setup_widget_class
        self.is_rendered = False
        self.is_loaded = False

    @property
    def is_ready(self):
        return self.is_loaded and self.is_active and self.selected

    def activate(self):
        self.is_active = True

    def deactivate(self):
        self.is_active = False

    def get_setup_widget(self) -> QEAppComputationalResourcesWidget:
        if not self.is_loaded:
            self._setup_widget = self.setup_widget_class(
                description=self.description,
                default_calc_job_plugin=self.default_calc_job_plugin,
            )
            self.is_loaded = True
        return self._setup_widget


CodesDict = dict[str, CodeModel]


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

    installing_qe = tl.Bool(False)
    installing_sssp = tl.Bool(False)
    qe_installed = tl.Bool(allow_none=True)
    sssp_installed = tl.Bool(allow_none=True)

    codes = tl.Dict(
        key_trait=tl.Unicode(),  # plugin identifier
        value_trait=tl.Dict(  # plugin codes
            key_trait=tl.Unicode(),  # code name
            value_trait=tl.Instance(CodeModel),  # code metadata
        ),
        default_value={},
    )

    def get_properties(self) -> list[str]:
        return self.input_parameters.get("workchain", {}).get("properties", [])

    def get_model_state(self) -> dict[str, dict[str, dict]]:
        return {
            "codes": self.get_selected_codes(),
        }

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

    # TODO identifier not the clearest name
    def add_code(self, identifier, name, code):
        if identifier not in self.codes:
            self.codes[identifier] = {}  # type: ignore
        self.codes[identifier][name] = code  # type: ignore

    def get_code(self, identifier, name) -> CodeModel | None:
        if identifier in self.codes and name in self.codes[identifier]:  # type: ignore
            return self.codes[identifier][name]  # type: ignore

    def get_codes(self, flat=False) -> t.Iterator[tuple[str, CodesDict | CodeModel]]:
        if flat:
            for codes in self.codes.values():
                yield from codes.items()
        else:
            yield from self.codes.items()

    def get_selected_codes(self) -> dict[str, dict]:
        return {
            name: code.parameters
            for name, code in self.get_codes(flat=True)
            if code.is_ready
        }  # type: ignore

    def set_selected_codes(self, code_data=DEFAULT["codes"]):
        def get_code_uuid(code):
            if code is not None:
                try:
                    return orm.load_code(code).uuid
                except NotExistent:
                    return None

        with self.hold_trait_notifications():
            for name, code in self.get_codes(flat=True):
                if name in code_data:
                    parameters = code_data[name]
                    code_uuid = get_code_uuid(parameters["code"])
                    if code_uuid in [opt[1] for opt in code.options]:
                        parameters["code"] = code_uuid
                        code.parameters = parameters

    def reset(self):
        with self.hold_trait_notifications():
            self.input_structure = None  # TODO why?
            self.process = None
            self.process_label = ""
            self.process_description = ""
            self.submission_blocker_messages = ""
