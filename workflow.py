"""Widgets for the submission of bands work chains.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
from pprint import pformat

import ipywidgets as ipw
import traitlets
from aiida.common import NotExistent
from aiida.engine import ProcessState
from aiida.engine import submit
from aiida.orm import ProcessNode
from aiida.orm import load_code
from aiida.plugins import DataFactory
from aiidalab_widgets_base import CodeDropdown
from aiidalab_widgets_base import ProcessMonitor
from aiidalab_widgets_base import ProcessNodesTreeWidget
from aiidalab_widgets_base import WizardAppWidgetStep

from pseudos import PseudoFamilySelector
from widgets import NodeViewWidget
from widgets import ResourceSelectionWidget

from apps.quantumespresso.qe_workflow import QeAppWorkChain
from aiida_quantumespresso.common.types import SpinType, ElectronicType, RelaxType


StructureData = DataFactory("structure")


def update_resources(builder, resources):
    for k, v in builder.items():
        if isinstance(v, dict):
            if k == "resources":
                builder["resources"].update(resources)
            else:
                update_resources(v, resources)


class WorkChainConfig(ipw.VBox):
    def __init__(self, **kwargs):

        self.run_geo_opt = ipw.Checkbox(description="Optimize geometry", value=True)
        self.geo_opt_type = ipw.Dropdown(
            description="Method:", options=["POSITIONS", "POSITIONS_CELL"]
        )
        ipw.dlink(
            (self.run_geo_opt, "value"),
            (self.geo_opt_type, "disabled"),
            transform=lambda v: not v,
        )

        self.run_bands = ipw.Checkbox(description="Calculate band structure")
        self.run_pdos = ipw.Checkbox(description="Calculate density of states")

        # Simulation protocol.
        self.simulation_protocol = ipw.Dropdown(
            options=["fast", "moderate", "precise"], value="moderate"
        )

        super().__init__(
            children=[
                ipw.HTML("Select which steps to run as part of the work chain:"),
                ipw.HBox(children=[self.run_geo_opt, self.geo_opt_type]),
                ipw.HBox(children=[self.run_bands]),
                self.run_pdos,
                ipw.HBox(
                    children=[
                        ipw.HTML("<b>Protocol</b>", layout=ipw.Layout(flex="1 1 auto")),
                        self.simulation_protocol,
                    ],
                    layout=ipw.Layout(vertical_align="bottom"),
                ),
                ipw.HTML(
                    """The "moderate" protocol represents a balanced trade-off between accuracy and speed.
            Choose the "fast" protocol for a faster calculation with less precision and the "precise" protocol
            that provides more accuracy but will take longer."""
                ),
            ],
            layout=ipw.Layout(max_width="600px"),
            **kwargs
        )


class OptionsConfig(ipw.VBox):
    spin_type = traitlets.Instance(SpinType, allow_none=True)
    electronic_type = traitlets.Instance(ElectronicType, allow_none=True)
    kpoints_distance = traitlets.Float(allow_none=True)
    degauss = traitlets.Float(allow_none=True)

    def __init__(self, **kwargs):

        # Spin type.
        self._set_spin_automatically = ipw.Checkbox(
            description="Use default spin type.",
            indent=False,
            value=True,
        )
        self._spin_type = ipw.Dropdown(
            options=["NONE", "COLLINEAR", "NON_COLLINEAR"],
            value="NONE",
            description="Spin Type:",
        )
        ipw.dlink(
            (self._set_spin_automatically, "value"), (self._spin_type, "disabled")
        )
        self._set_spin_automatically.observe(self.set_spin_type_trait, "value")
        self._spin_type.observe(self.set_spin_type_trait, "value")
        self.set_spin_type_trait()

        # Electronic type.
        self._set_el_type_automatically = ipw.Checkbox(
            description="Use default electronic type.",
            indent=False,
            value=True,
        )
        self._electronic_type = ipw.Dropdown(
            options=["METAL", "INSULATOR"],
            value="METAL",
            description="Electronic Type:",
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self._set_el_type_automatically, "value"),
            (self._electronic_type, "disabled"),
        )
        self._set_el_type_automatically.observe(self.set_electronic_type_trait, "value")
        self._electronic_type.observe(self.set_electronic_type_trait, "value")
        self.set_electronic_type_trait()

        # K-points distance.
        self._set_kpoints_distance_automatically = ipw.Checkbox(
            description="Use default k-points distance.",
            indent=False,
            value=True,
        )
        self._kpoints_distance = ipw.FloatText(
            value=0.4,
            step=0.1,
            description="K-points distance:",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self._set_kpoints_distance_automatically, "value"),
            (self._kpoints_distance, "disabled"),
        )
        self._kpoints_distance.observe(self.set_kpoints_distance_trait, "value")
        self._set_kpoints_distance_automatically.observe(
            self.set_kpoints_distance_trait, "value"
        )
        self.set_kpoints_distance_trait()

        # Modify degauss.
        # degauss = ipw.FloatText(
        #    value=0.01,
        #    step=0.01,
        #    description="Gaussian spreading",
        #    disabled=False,
        #    style={"description_width": "initial"},
        # )

        super().__init__(
            children=[
                ipw.HBox([self._set_spin_automatically, self._spin_type]),
                ipw.HBox([self._set_el_type_automatically, self._electronic_type]),
                ipw.HBox(
                    [self._set_kpoints_distance_automatically, self._kpoints_distance]
                ),
            ],
            layout=ipw.Layout(max_width="600px"),
            **kwargs
        )

    _DEFAULT_SPIN_TYPE = "NONE"

    def set_spin_type_trait(self, _=None):
        self.spin_type = SpinType[
            self._DEFAULT_SPIN_TYPE
            if self._set_spin_automatically.value
            else self._spin_type.value.upper()
        ]

    @traitlets.observe("spin_type")
    def _observe_spin_type(self, change):
        with self.hold_trait_notifications():
            self._spin_type.value = change["new"].value.upper()
            self._set_spin_automatically.value = (
                change["new"] is self._DEFAULT_SPIN_TYPE
            )

    _DEFAULT_ELECTRONIC_TYPE = "METAL"

    def set_electronic_type_trait(self, _=None):
        self.electronic_type = ElectronicType[
            self._DEFAULT_ELECTRONIC_TYPE
            if self._set_el_type_automatically.value
            else self._electronic_type.value
        ]

    @traitlets.observe("electronic_type")
    def _observe_electronic_type(self, change):
        with self.hold_trait_notifications():
            self._electronic_type.value = change["new"].value.upper()
            self._set_el_type_automatically.value = (
                change["new"] is self._DEFAULT_ELECTRONIC_TYPE
            )

    _DEFAULT_KPOINTS_DISTANCE = None

    def set_kpoints_distance_trait(self, _=None):
        self.kpoints_distance = (
            self._DEFAULT_KPOINTS_DISTANCE
            if self._set_kpoints_distance_automatically.value
            else self._kpoints_distance.value
        )

    @traitlets.observe("kpoints_distance")
    def _observe_kpoints_distance(self, change):
        with self.hold_trait_notifications():
            self._kpoints_distance.value = (
                0.4 if change["new"] is None else change["new"].value.upper()
            )
            self._set_kpoints_distance_automatically.value = (
                change["new"] is self._DEFAULT_KPOINTS_DISTANCE
            )


class CodesConfig(ipw.VBox):
    def __init__(self, **kwargs):

        self.pw = CodeDropdown(
            input_plugin="quantumespresso.pw",
            description="PW code:",
            setup_code_params={
                "computer": "localhost",
                "description": "pw.x in AiiDAlab container.",
                "label": "pw",
                "input_plugin": "quantumespresso.pw",
                "remote_abs_path": "/usr/bin/pw.x",
            },
        )
        self.dos = CodeDropdown(
            input_plugin="quantumespresso.dos",
            description="DOS code",
            setup_code_params={
                "computer": "localhost",
                "description": "dos.x in AiiDAlab container.",
                "label": "dos",
                "input_plugin": "quantumespresso.dos",
                "remote_abs_path": "/usr/bin/dos.x",
            },
        )

        self.projwfc = CodeDropdown(
            input_plugin="quantumespresso.projwfc",
            description="PROJWFC code",
            setup_code_params={
                "computer": "localhost",
                "description": "projwfc.x in AiiDAlab container.",
                "label": "projwfc",
                "input_plugin": "quantumespresso.projwfc",
                "remote_abs_path": "/usr/bin/projwfc.x",
            },
        )
        super().__init__(
            children=[
                self.pw,
                self.dos,
                self.projwfc,
            ],
            **kwargs
        )


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    input_structure = traitlets.Instance(StructureData, allow_none=True)
    process = traitlets.Instance(ProcessNode, allow_none=True)
    disabled = traitlets.Bool()
    builder_parameters = traitlets.Dict()
    expert_mode = traitlets.Bool()

    def __init__(self, **kwargs):
        self.workchain_config = WorkChainConfig()
        self.resources_config = ResourceSelectionWidget()
        self.options_config = OptionsConfig()
        self.pseudo_family_selector = PseudoFamilySelector()
        self.codes_selector = CodesConfig()

        self.codes_selector.pw.observe(self._update_state, "selected_code")
        self.codes_selector.dos.observe(self._update_state, "selected_code")
        self.codes_selector.projwfc.observe(self._update_state, "selected_code")
        self.workchain_config.run_pdos.observe(self._update_state, "value")

        self._setup_builder_parameters_update()

        self.tab = ipw.Tab(
            children=[
                self.workchain_config,
                self.resources_config,
                self.codes_selector,
            ],
            layout=ipw.Layout(min_height="250px"),
        )

        self.tab.set_title(0, "Workchain")
        self.tab.set_title(1, "Compute resources")
        self.tab.set_title(2, "Select codes")

        self.submit_button = ipw.Button(
            description="Submit",
            tooltip="Submit the calculation with the selected parameters.",
            icon="play",
            button_style="success",
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
            disabled=True,
        )

        self.submit_button.on_click(self._on_submit_button_clicked)

        self.expert_mode_control = ipw.Checkbox(description="Expert mode")
        ipw.link((self, "expert_mode"), (self.expert_mode_control, "value"))

        self._update_builder_parameters()

        self.builder_parameters_view = ipw.HTML(layout=ipw.Layout(width="auto"))
        ipw.dlink(
            (self, "builder_parameters"),
            (self.builder_parameters_view, "value"),
            transform=lambda p: '<pre style="line-height: 100%">'
            + pformat(p, indent=2, width=200)
            + "</pre>",
        )

        super().__init__(
            children=[
                self.tab,
                ipw.HBox([self.submit_button, self.expert_mode_control]),
            ]
        )

    @traitlets.observe("expert_mode")
    def _observe_expert_mode(self, change):
        if change["new"]:
            self.tab.set_title(0, "Workchain")
            self.tab.set_title(1, "Compute resources")
            self.tab.set_title(2, "Advanced options")
            self.tab.set_title(3, "Pseudo potentials")
            self.tab.set_title(4, "Select codes")
            self.tab.set_title(5, "Parameters")
            self.tab.children = [
                self.workchain_config,
                self.resources_config,
                self.options_config,
                self.pseudo_family_selector,
                self.codes_selector,
                self.builder_parameters_view,
            ]
        else:
            self.tab.set_title(0, "Workchain")
            self.tab.set_title(1, "Compute resources")
            self.tab.set_title(2, "Select codes")
            self.tab.children = [
                self.workchain_config,
                self.resources_config,
                self.codes_selector,
            ]

    def _get_state(self):

        # Process is already running.
        if self.process is not None:
            return self.State.SUCCESS

        # Input structure not specified.
        if self.input_structure is None:
            return self.State.INIT

        # Pseudo family is not installed.
        if not self.pseudo_family_selector.installed:
            return self.State.READY

        # PW code not selected.
        if self.codes_selector.pw.selected_code is None:
            return self.State.READY

        # PDOS run requested, but codes are not specified.
        if self.workchain_config.run_pdos.value:
            if self.codes_selector.dos.selected_code is None:
                return self.State.READY
            if self.codes_selector.projwfc.selected_code is None:
                return self.State.READY

        return self.State.CONFIGURED

    def _update_state(self, _=None):
        self.state = self._get_state()

    @traitlets.observe("state")
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            self.disabled = change["new"] not in (
                self.State.READY,
                self.State.CONFIGURED,
            )
            self.submit_button.disabled = change["new"] != self.State.CONFIGURED

    @traitlets.observe("input_structure")
    def _observe_input_structure(self, change):
        self.set_trait("builder_parameters", self._default_builder_parameters())
        self._update_state()

    @traitlets.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            if process_node is not None:
                self.input_structure = process_node.inputs.structure
                builder_parameters = process_node.get_extra("builder_parameters", None)
                if builder_parameters is not None:
                    self.set_trait("builder_parameters", builder_parameters)
            self._update_state()

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.submit()

    def _setup_builder_parameters_update(self):
        update = self._update_builder_parameters  # alias for code conciseness
        # Codes
        self.codes_selector.dos.observe(update, ["selected_code"])
        self.codes_selector.projwfc.observe(update, ["selected_code"])
        self.codes_selector.pw.observe(update, ["selected_code"])
        # Protocol and additional parameters
        self.options_config.observe(update, ["electronic_type"])
        self.options_config.observe(update, ["kpoints_distance"])
        self.options_config.observe(update, ["spin_type"])
        self.pseudo_family_selector.observe(update, ["value"])
        self.workchain_config.geo_opt_type.observe(update, ["value"])
        self.workchain_config.run_geo_opt.observe(update, ["value"])
        self.workchain_config.simulation_protocol.observe(update, ["value"])
        # "Extra" parameters
        self.workchain_config.run_bands.observe(update, ["value"])
        self.workchain_config.run_pdos.observe(update, ["value"])

    @staticmethod
    def _serialize_builder_parameters(parameters):
        parameters = parameters.copy()  # create copy to not modify original dict

        # Codes
        def _get_uuid(code):
            return None if code is None else str(code.uuid)

        parameters["dos_code"] = _get_uuid(parameters["dos_code"])
        parameters["projwfc_code"] = _get_uuid(parameters["projwfc_code"])
        parameters["pw_code"] = _get_uuid(parameters["pw_code"])

        # Protocol
        parameters["electronic_type"] = parameters["electronic_type"].value
        parameters["relax_type"] = parameters["relax_type"].value
        parameters["spin_type"] = parameters["spin_type"].value
        return parameters

    @staticmethod
    def _deserialize_builder_parameters(parameters):
        parameters = parameters.copy()  # create copy to not modify original dict

        # Codes
        def _load_code(code):
            if code is not None:
                try:
                    return load_code(code)
                except NotExistent as error:
                    print("error", error)
                    return None

        parameters["dos_code"] = _load_code(parameters["dos_code"])
        parameters["projwfc_code"] = _load_code(parameters["projwfc_code"])
        parameters["pw_code"] = _load_code(parameters["pw_code"])

        # Protocol
        parameters["electronic_type"] = ElectronicType(parameters["electronic_type"])
        parameters["relax_type"] = RelaxType(parameters["relax_type"])
        parameters["spin_type"] = SpinType(parameters["spin_type"])
        return parameters

    def _update_builder_parameters(self, _=None):
        self.set_trait(
            "builder_parameters",
            self._serialize_builder_parameters(
                dict(
                    # Codes
                    dos_code=self.codes_selector.dos.selected_code,
                    projwfc_code=self.codes_selector.projwfc.selected_code,
                    pw_code=self.codes_selector.pw.selected_code,
                    # Protocol and additional parameters
                    electronic_type=self.options_config.electronic_type,
                    kpoints_distance_override=self.options_config.kpoints_distance,
                    protocol=self.workchain_config.simulation_protocol.value,
                    pseudo_family=self.pseudo_family_selector.value,
                    relax_type=RelaxType[self.workchain_config.geo_opt_type.value]
                    if self.workchain_config.run_geo_opt.value
                    else RelaxType["NONE"],
                    spin_type=self.options_config.spin_type,
                    # "extra" parameters
                    run_bands=self.workchain_config.run_bands.value,
                    run_pdos=self.workchain_config.run_pdos.value,
                )
            ),
        )

    @traitlets.observe("builder_parameters")
    def _observe_builder_parameters(self, change):
        bp = self._deserialize_builder_parameters(change["new"])

        with self.hold_trait_notifications():
            # Codes
            self.codes_selector.dos.selected_code = bp.get("dos_code")
            self.codes_selector.projwfc.selected_code = bp.get("projwfc_code")
            self.codes_selector.pw.selected_code = bp.get("pw_code")
            # Protocol and additional parameters
            self.options_config.electronic_type = bp.get(
                "electronic_type", ElectronicType["METAL"]
            )
            self.options_config.kpoints_distance = bp.get("kpoints_distance")
            self.workchain_config.simulation_protocol.value = bp["protocol"]
            self.pseudo_family_selector.value = bp["pseudo_family"]

            relax_type = bp.get("relax_type", RelaxType["NONE"])
            self.workchain_config.geo_opt_type.value = (
                "POSITIONS"
                if relax_type is RelaxType["NONE"]
                else relax_type.value.upper()
            )
            self.workchain_config.run_geo_opt.value = (
                relax_type is not RelaxType["NONE"]
            )

            self.options_config.spin_type = bp.get("spin_type", SpinType["NONE"])
            # "extra" parameters
            self.workchain_config.run_bands.value = bp["run_bands"]
            self.workchain_config.run_pdos.value = bp["run_pdos"]

    def submit(self, _=None):

        assert self.input_structure is not None

        builder_parameters = self.builder_parameters.copy()
        builder_parameters_for_extras = builder_parameters.copy()

        run_bands = builder_parameters.pop("run_bands")
        run_pdos = builder_parameters.pop("run_pdos")

        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            **self._deserialize_builder_parameters(builder_parameters)
        )

        if not run_bands:
            builder.pop("bands")
        if not run_pdos:
            builder.pop("pdos")

        resources = {
            "num_machines": self.resources_config.number_of_nodes.value,
            "num_mpiprocs_per_machine": self.resources_config.cpus_per_node.value,
        }
        update_resources(builder, resources)

        self.process = submit(builder)
        self.process.set_extra("builder_parameters", builder_parameters_for_extras)

    def reset(self):
        with self.hold_trait_notifications():
            self.process = None
            self.input_structure = None
            self.builder_parameters = self._default_builder_parameters()

    @traitlets.default("builder_parameters")
    def _default_builder_parameters(self):
        return {
            # Codes
            "dos_code": "dos@localhost",
            "projwfc_code": "projwfc@localhost",
            "pw_code": "pw@localhost",
            # Protocol and additional parameters
            "electronic_type": "metal",
            "kpoints_distance_override": None,
            "protocol": "moderate",
            "pseudo_family": None,
            "relax_type": "positions",
            "spin_type": "none",
            # Extra parameters
            "run_bands": False,
            "run_pdos": False,
        }


class ViewQeAppWorkChainStatusAndResultsStep(ipw.VBox, WizardAppWidgetStep):

    process = traitlets.Instance(ProcessNode, allow_none=True)

    def __init__(self, **kwargs):
        self.process_tree = ProcessNodesTreeWidget()
        ipw.dlink((self, "process"), (self.process_tree, "process"))

        self.node_view = NodeViewWidget(
            layout={"width": "auto", "height": "auto", "border": "1px solid black"}
        )
        ipw.dlink(
            (self.process_tree, "selected_nodes"),
            (self.node_view, "node"),
            transform=lambda nodes: nodes[0] if nodes else None,
        )
        self.process_status = ipw.VBox(children=[self.process_tree, self.node_view])

        # Setup process monitor
        self.process_monitor = ProcessMonitor(
            timeout=0.2,
            callbacks=[
                self.process_tree.update,
                self._update_state,
            ],
        )
        ipw.dlink((self, "process"), (self.process_monitor, "process"))

        super().__init__([self.process_status], **kwargs)

    def can_reset(self):
        "Do not allow reset while process is running."
        return self.state is not self.State.ACTIVE

    def reset(self):
        self.process = None

    def _update_state(self):
        if self.process is None:
            self.state = self.State.INIT
        else:
            process_state = self.process.process_state
            if process_state in (
                ProcessState.CREATED,
                ProcessState.RUNNING,
                ProcessState.WAITING,
            ):
                self.state = self.State.ACTIVE
            elif process_state in (ProcessState.EXCEPTED, ProcessState.KILLED):
                self.state = self.State.FAIL
            elif process_state is ProcessState.FINISHED:
                self.state = self.State.SUCCESS

    @traitlets.observe("process")
    def _observe_process(self, change):
        self._update_state()
