from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl
from IPython.display import Javascript, display

from aiida import orm
from aiida.engine import ProcessBuilderNamespace, submit
from aiida.orm.utils.serialize import serialize
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.common.mixins import HasInputStructure, HasModels, HasProcess
from aiidalab_qe.common.panel import PluginResourceSettingsModel, ResourceSettingsModel
from aiidalab_qe.common.wizard import ConfirmableDependentWizardStepModel, State
from aiidalab_qe.utils import shallow_copy_nested_dict
from aiidalab_qe.workflows import QeAppWorkChain

DEFAULT: dict = DEFAULT_PARAMETERS  # type: ignore


class SubmissionStepModel(
    ConfirmableDependentWizardStepModel,
    HasModels[ResourceSettingsModel],
    HasInputStructure,
    HasProcess,
):
    identifier = "submission"

    input_parameters = tl.Dict()

    process_label = tl.Unicode("")
    process_description = tl.Unicode("")

    warning_messages = tl.Unicode("")

    installing_qe = tl.Bool(False)
    qe_installed = tl.Bool(allow_none=True)

    plugin_overrides = tl.List(tl.Unicode())

    fetched_resources = tl.Bool(False)

    _dependencies = [
        "structure_uuid",
        "input_parameters",
    ]

    def __init__(self, *args, **kwargs):
        self._default_models = {
            "global",
        }

        super().__init__(*args, **kwargs)

        self.confirmation_exceptions += [
            "warning_messages",
            "installing_qe",
            "qe_installed",
            "fetched_resources",
        ]

        self.default_user_email = orm.User.collection.get_default().email

    def confirm(self):
        super().confirm()
        if not self.has_process:
            self._submit()

    def update(self):
        self.update_process_label()
        self.update_plugin_inclusion()
        self.update_plugin_overrides()
        self.update_blockers()
        for _, model in self.get_models():
            model.update()

    def update_process_label(self):
        if not self.has_structure or not self.input_parameters:
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
        if not self.input_parameters:
            return
        properties = self._get_properties()
        for identifier, model in self.get_models():
            if identifier in self._default_models:
                continue
            model.include = identifier in properties

    def update_plugin_overrides(self):
        if not self.input_parameters:
            return
        self.plugin_overrides = [
            identifier
            for identifier, model in self.get_models()
            if isinstance(model, PluginResourceSettingsModel)
            and model.include
            and model.override
        ]

    def update_warnings(self):
        warning_messages = self._check_warnings()
        for _, model in self.get_models():
            warning_messages += model.warning_messages
        self.warning_messages = warning_messages

    def update_process_metadata(self):
        if not self.has_process:
            return
        self.process_label = self.process.label
        self.process_description = self.process.description
        self.locked = True

    def get_model_state(self) -> dict:
        return {
            identifier: model.get_model_state()
            for identifier, model in self.get_models()
            if model.include
        }

    def set_model_state(self, state: dict):
        for identifier, model in self.get_models():
            if state.get(identifier):
                model.include = True
                model.set_model_state(state[identifier])

    def update_state(self):
        if self.confirmed:
            self.state = State.SUCCESS
        elif self.is_previous_step_successful and (
            self.has_structure and self.input_parameters
        ):
            self.state = State.CONFIGURED
        else:
            self.state = State.INIT

    def reset(self):
        with self.hold_trait_notifications():
            self.input_parameters = {}
            self.process_uuid = None
            for identifier, model in self.get_models():
                if identifier not in self._default_models:
                    model.include = False

    def _submit(self):
        parameters = shallow_copy_nested_dict(self.input_parameters)
        parameters |= {"codes": self.get_model_state()}
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
            self.process_uuid = process_node.uuid

            pk = process_node.pk
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
            # Use shallow copy to avoid parameter mutation by the workflow
            parameters=shallow_copy_nested_dict(parameters),
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

    def _check_blockers(self):
        if not self.has_structure:
            yield "No selected input structure"

        if not self.input_parameters:
            yield "No selected input parameters"

        if self.installing_qe:
            yield "Installing Quantum ESPRESSO codes..."

        if not self.qe_installed:
            yield "Quantum ESPRESSO is not yet installed"

    def _check_warnings(self):
        """Check for any warnings that should be displayed to the user."""
        return ""
