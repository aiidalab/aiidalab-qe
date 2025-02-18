from __future__ import annotations

from copy import deepcopy

import ipywidgets as ipw
import traitlets as tl
from IPython.display import Javascript, display

from aiida import orm
from aiida.engine import ProcessBuilderNamespace, submit
from aiida.orm.utils.serialize import serialize
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.mixins import Confirmable, HasInputStructure, HasModels
from aiidalab_qe.common.panel import PluginResourceSettingsModel, ResourceSettingsModel
from aiidalab_qe.common.widgets import QeWizardStepModel
from aiidalab_qe.workflows import QeAppWorkChain

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmissionStepModel(
    QeWizardStepModel,
    HasModels[ResourceSettingsModel],
    HasInputStructure,
    Confirmable,
):
    identifier = "submission"

    input_parameters = tl.Dict()

    process_node = tl.Instance(orm.WorkChainNode, allow_none=True)
    process_label = tl.Unicode("")
    process_description = tl.Unicode("")

    internal_submission_blockers = tl.List(tl.Unicode())
    external_submission_blockers = tl.List(tl.Unicode())
    submission_blocker_messages = tl.Unicode("")
    submission_warning_messages = tl.Unicode("")

    plugin_overrides = tl.List(tl.Unicode())

    confirmation_exceptions = [
        "confirmed",
        "internal_submission_blockers",
        "external_submission_blockers",
        "submission_blocker_messages",
        "submission_warning_messages",
    ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self._default_models = {
            "global",
        }

        self._ALERT_MESSAGE = """
            <div class="alert alert-{alert_class} alert-dismissible">
                <a href="#" class="close" data-dismiss="alert" aria-label="close">
                    &times;
                </a>
                <strong>{message}</strong>
            </div>
        """

    @property
    def is_blocked(self):
        return any(
            [
                *self.internal_submission_blockers,
                *self.external_submission_blockers,
            ]
        )

    def confirm(self):
        super().confirm()
        if not self.process_node:
            self._submit()

    def update(self):
        for _, model in self.get_models():
            model.update()

    def update_process_label(self):
        if not self.input_structure:
            self.process_label = ""
            return
        structure_label = (
            self.input_structure.label
            if len(self.input_structure.label) > 0
            else self.input_structure.get_formula()
        )
        workchain_data = self.input_parameters.get(
            "workchain",
            {"properties": []},
        )

        soc_parameters = (
            self.input_parameters.get("advanced", {})
            .get("pw", {})
            .get("parameters", {})
            .get("SYSTEM", {})
            .get("lspinorb", False)
        )

        soc_info = "spin-orbit coupling" if soc_parameters else ""

        properties = [p for p in workchain_data["properties"] if p != "relax"]
        relax_type = workchain_data.get("relax_type", "none")
        relax_info = "unrelaxed"
        if relax_type != "none":
            relax_info = (
                "relax: atoms+cell" if "cell" in relax_type else "relax: atoms only"
            )
        protocol_and_magnetic_info = f"{workchain_data['protocol']} protocol"
        if workchain_data["spin_type"] != "none":
            protocol_and_magnetic_info += ", magnetic"
        properties_info = f"→ {', '.join(properties)}" if properties else ""

        label_details = [
            relax_info,
            protocol_and_magnetic_info,
            soc_info,
        ]
        filtered_label_details = [detail for detail in label_details if detail]
        label = f"{structure_label} [{', '.join(filtered_label_details)}] {properties_info}".strip()

        self.process_label = label

    def update_plugin_inclusion(self):
        properties = self._get_properties()
        for identifier, model in self.get_models():
            if identifier in self._default_models:
                continue
            model.include = identifier in properties

    def update_plugin_overrides(self):
        self.plugin_overrides = [
            identifier
            for identifier, model in self.get_models()
            if isinstance(model, PluginResourceSettingsModel)
            and model.include
            and model.override
        ]

    def update_submission_blockers(self):
        submission_blockers = []
        for _, model in self.get_models():
            submission_blockers += model.submission_blockers
        self.internal_submission_blockers = submission_blockers

    def update_submission_warnings(self):
        self.submission_warning_messages = "".join(
            model.submission_warning_messages for _, model in self.get_models()
        )  # type: ignore

    def update_submission_blocker_message(self):
        blockers = self.internal_submission_blockers + self.external_submission_blockers
        if any(blockers):
            formatted = "\n".join(f"<li>{item}</li>" for item in blockers)
            self.submission_blocker_messages = f"""
                <div class="alert alert-danger">
                    <b>The submission is blocked due to the following reason(s):</b>
                    <ul>
                        {formatted}
                    </ul>
                </div>
            """
        else:
            self.submission_blocker_messages = ""

    def get_model_state(self) -> dict[str, dict[str, dict]]:
        parameters: dict = deepcopy(self.input_parameters)  # type: ignore
        parameters["codes"] = {
            identifier: model.get_model_state()
            for identifier, model in self.get_models()
            if model.include
        }
        return parameters

    def set_model_state(self, parameters):
        codes: dict = parameters.get("codes", {})

        if "resources" in parameters:
            resources = parameters["resources"]
            codes |= {key: {"code": value} for key, value in codes.items()}
            codes["pw"]["nodes"] = resources["num_machines"]
            codes["pw"]["cpus"] = resources["num_mpiprocs_per_machine"]
            codes["pw"]["parallelization"] = {"npool": resources["npools"]}

        workchain_parameters: dict = parameters.get("workchain", {})
        properties = set(workchain_parameters.get("properties", []))
        included = self._default_models | properties
        for identifier, model in self.get_models():
            model.include = identifier in included
            if codes.get(identifier):
                model.set_model_state(codes[identifier])
                model.loaded_from_process = True

        if self.process_node:
            self.process_label = self.process_node.label
            self.process_description = self.process_node.description
            self.loaded_from_process = True

    def get_selected_codes(self) -> dict[str, dict]:
        return {
            identifier: code_model.get_model_state()
            for identifier, code_model in self.get_model("global").get_models()
            if code_model.is_ready
        }

    def reset(self):
        with self.hold_trait_notifications():
            self.input_structure = None
            self.input_parameters = {}
            self.process_node = None
            for identifier, model in self.get_models():
                if identifier not in self._default_models:
                    model.include = False

    def _submit(self):
        parameters = self.get_model_state()
        builder = self._create_builder(parameters)

        with self.hold_trait_notifications():
            process_node = submit(builder)

            process_node.label = self.process_label
            process_node.description = self.process_description
            # since AiiDA data node may exist in the ui_parameters,
            # we serialize it to yaml
            process_node.base.extras.set("ui_parameters", serialize(parameters))
            # store the workchain name in extras, this will help to filter the workchain in the future
            process_node.base.extras.set("workchain", parameters["workchain"])  # type: ignore
            process_node.base.extras.set(
                "structure",
                self.input_structure.get_formula(),
            )
            self.process_node = process_node

            self._update_url()

    def _update_url(self):
        pk = self.process_node.pk
        display(Javascript(f"window.history.pushState(null, '', '?pk={pk}');"))

    def _link_model(self, model: ResourceSettingsModel):
        for dependency in model.dependencies:
            dependency_parts = dependency.split(".")
            if len(dependency_parts) == 1:  # from parent, e.g. input_structure
                target_model = self
                trait = dependency
            else:  # from sibling, e.g. workchain.protocol
                sibling, trait = dependency_parts
                target_model = self.get_model(sibling)
            ipw.dlink(
                (target_model, trait),
                (model, trait),
            )

    def _get_properties(self) -> list[str]:
        return self.input_parameters.get("workchain", {}).get("properties", [])

    def _create_builder(self, parameters) -> ProcessBuilderNamespace:
        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            parameters=deepcopy(parameters),  # TODO why deepcopy again?
        )

        codes = parameters["codes"]["global"]["codes"]
        if "relax" in builder:
            builder.relax.base.pw.metadata.options.resources = {
                "num_machines": codes.get("quantumespresso__pw")["nodes"],
                "num_mpiprocs_per_machine": codes.get("quantumespresso__pw")[
                    "ntasks_per_node"
                ],
                "num_cores_per_mpiproc": codes.get("quantumespresso__pw")[
                    "cpus_per_task"
                ],
            }
            mws = codes.get("quantumespresso__pw")["max_wallclock_seconds"]
            builder.relax.base.pw.metadata.options["max_wallclock_seconds"] = mws
            parallelization = codes["quantumespresso__pw"]["parallelization"]
            builder.relax.base.pw.parallelization = orm.Dict(dict=parallelization)

        return builder
