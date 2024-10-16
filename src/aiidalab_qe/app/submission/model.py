from __future__ import annotations

import typing as t
from copy import deepcopy

import traitlets as tl

from aiida import orm
from aiida.common import NotExistent
from aiida.engine import ProcessBuilderNamespace
from aiida.engine import submit as aiida_submit
from aiida.orm.utils.serialize import serialize
from aiida_quantumespresso.data.hubbard_structure import HubbardStructureData
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.widgets import QEAppComputationalResourcesWidget
from aiidalab_qe.workflows import QeAppWorkChain

from .code import CodeModel, CodesDict

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

    internal_submission_blockers = tl.List(tl.Unicode())
    external_submission_blockers = tl.List(tl.Unicode())

    code_widgets: dict[str, QEAppComputationalResourcesWidget] = {}

    @property
    def is_blocked(self):
        return any(
            [
                *self.internal_submission_blockers,
                *self.external_submission_blockers,
            ]
        )

    def submit(self):
        parameters = self._get_submission_parameters()
        builder = self._create_builder(parameters)

        with self.hold_trait_notifications():
            process = aiida_submit(builder)

            process.label = self.process_label
            process.description = self.process_description
            # since AiiDA data node may exist in the ui_parameters,
            # we serialize it to yaml
            process.base.extras.set("ui_parameters", serialize(parameters))
            # store the workchain name in extras, this will help to filter the workchain in the future
            process.base.extras.set("workchain", parameters["workchain"])  # type: ignore
            process.base.extras.set(
                "structure",
                self.input_structure.get_formula(),
            )
            self.process = process

    def update_active_codes(self):
        for name, code in self.get_codes(flat=True):
            if name != "pw":
                code.deactivate()
        properties = self.get_properties()
        for identifier, codes in self.get_codes():
            if identifier in properties:
                for code in codes.values():
                    code.activate()

    def update_process_label(self):
        if not self.input_structure:
            self.process_label = ""
            return
        formula = self.input_structure.get_formula()
        workchain_data = self.input_parameters.get(
            "workchain",
            {"properties": []},
        )
        properties = [p for p in workchain_data["properties"] if p != "relax"]
        relax_type = workchain_data.get("relax_type", "none")
        if relax_type != "none":
            relax_info = "structure is relaxed"
        else:
            relax_info = "structure is not relaxed"
        if not properties:
            properties_info = ""
        else:
            properties_info = f", properties on {', '.join(properties)}"

        label = f"{formula} {relax_info} {properties_info}"
        self.process_label = label

    def update_submission_blockers(self):
        self.internal_submission_blockers = list(self._check_submission_blockers())

    def update_submission_blocker_message(self):
        blockers = self.internal_submission_blockers + self.external_submission_blockers
        if any(blockers):
            fmt_list = "\n".join(f"<li>{item}</li>" for item in sorted(blockers))
            self.submission_blocker_messages = f"""
                <div class="alert alert-info">
                    <b>The submission is blocked due to the following reason(s):</b>
                    <ul>
                        {fmt_list}
                    </ul>
                </div>
            """
        else:
            self.submission_blocker_messages = ""

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

    def add_code(self, identifier, name, code):
        code.name = name
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
            self.input_structure = None
            self.input_parameters = {}
            self.process = None

    def _create_builder(self, parameters) -> ProcessBuilderNamespace:
        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            parameters=deepcopy(parameters),  # TODO why deepcopy again?
        )

        codes = parameters["codes"]

        builder.relax.base.pw.metadata.options.resources = {
            "num_machines": codes.get("pw")["nodes"],
            "num_mpiprocs_per_machine": codes.get("pw")["ntasks_per_node"],
            "num_cores_per_mpiproc": codes.get("pw")["cpus_per_task"],
        }
        mws = codes.get("pw")["max_wallclock_seconds"]
        builder.relax.base.pw.metadata.options["max_wallclock_seconds"] = mws
        parallelization = codes["pw"]["parallelization"]
        builder.relax.base.pw.parallelization = orm.Dict(dict=parallelization)

        return builder

    def _get_submission_parameters(self) -> dict:
        submission_parameters = self.get_model_state()
        for name, code_widget in self.code_widgets.items():
            if name in submission_parameters["codes"]:
                for key, value in code_widget.parameters.items():
                    if key != "code":
                        submission_parameters["codes"][name][key] = value
        parameters = deepcopy(self.input_parameters)
        parameters.update(submission_parameters)
        return parameters  # type: ignore

    def _check_submission_blockers(self):
        # Do not submit while any of the background setup processes are running.
        if self.installing_qe or self.installing_sssp:
            yield "Background setup processes must finish."

        # SSSP library not installed
        if not self.sssp_installed:
            yield "The SSSP library is not installed."

        # No pw code selected (this is ignored while the setup process is running).
        pw_code = self.get_code(identifier="dft", name="pw")
        if pw_code and not pw_code.selected and not self.installing_qe:
            yield ("No pw code selected")

        # code related to the selected property is not installed
        properties = self.get_properties()
        message = "Calculating the {property} property requires code {code} to be set."
        for identifier, codes in self.get_codes():
            if identifier in properties:
                for code in codes.values():
                    if not code.is_ready:
                        yield message.format(property=identifier, code=code.description)

        # check if the QEAppComputationalResourcesWidget is used
        for name, code in self.get_codes(flat=True):
            # skip if the code is not displayed, convenient for the plugin developer
            if not code.is_ready:
                continue
            if not issubclass(
                code.setup_widget_class, QEAppComputationalResourcesWidget
            ):
                yield (
                    f"Error: hi, plugin developer, please use the QEAppComputationalResourcesWidget from aiidalab_qe.common.widgets for code {name}."
                )
