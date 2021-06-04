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
from aiida_quantumespresso.common.types import SpinType, ElectronicType, RelaxType
from aiidalab_qe_workchain import QeAppWorkChain
from aiidalab_widgets_base import CodeDropdown
from aiidalab_widgets_base import ProcessMonitor
from aiidalab_widgets_base import ProcessNodesTreeWidget
from aiidalab_widgets_base import WizardAppWidgetStep

from aiidalab_qe.pseudos import PseudoFamilySelector
from aiidalab_qe.widgets import NodeViewWidget
from aiidalab_qe.widgets import ResourceSelectionWidget


StructureData = DataFactory("structure")


def update_resources(builder, resources):
    for k, v in builder.items():
        if isinstance(v, dict):
            if k == "resources":
                builder["resources"].update(resources)
            else:
                update_resources(v, resources)


class WorkChainConfig(ipw.VBox):

    properties_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Properties</h4></div>"""
    )
    properties_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px">
        By default, the work chain will optimize the provided geometry. Uncheck
        the box if this is not desired.  The "POSITIONS" method will only
        optimize the atomic positions, "POSITIONS_CELL" will also optimize the
        unit cell of the structure.  The band structure work chain will
        automatically detect the default path in reciprocal space using the <a
        href="https://www.materialscloud.org/work/tools/seekpath target="_blank"> SeeK-path
        tool</a>.</div>"""
    )

    protocol_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Protocol</h4></div>"""
    )
    protocol_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 6px; padding-bottom: 0px">
        The "moderate" protocol represents a balanced trade-off between
        accuracy and speed. Choose the "fast" protocol for a faster calculation
        with less precision and the "precise" protocol that provides more
        accuracy but will take longer.</div>"""
    )

    def __init__(self, **kwargs):

        self.geo_opt_run = ipw.ToggleButton(value=True, description="Geometry", button_style='info')
        self.geo_opt_type = ipw.Dropdown(
            description="Method:",
            value="POSITIONS_CELL",
            options=["POSITIONS_CELL", "POSITIONS"],
        )
        ipw.dlink(
            (self.geo_opt_run, "value"),
            (self.geo_opt_type, "disabled"),
            transform=lambda v: not v,
        )
        self.bands_run = ipw.ToggleButton(indent=False, description="Band structure", button_style='info')
    
        self.run_pdos = ipw.ToggleButton( 
            description="Calculate density of states",
        )
        
        # Simulation protocol.
        self.simulation_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"], value="moderate", button_style='info'
        )
        super().__init__(
            children=[
                self.properties_title,
                ipw.HTML("Select which properties to calculate:"),
                ipw.HBox(
                    children=[self.geo_opt_run, self.geo_opt_type],
                ),
                self.bands_run,
                self.properties_help,
                self.protocol_title,
                ipw.HTML("Select the protocol:", layout=ipw.Layout(flex="1 1 auto")),
                self.simulation_protocol,
                self.protocol_help,
            ],
            **kwargs
        )


class OptionsConfig(ipw.VBox):
    spin_type = traitlets.Instance(SpinType, allow_none=True)
    electronic_type = traitlets.Instance(ElectronicType, allow_none=True)
    kpoints_distance = traitlets.Float(allow_none=True)
    degauss = traitlets.Float(allow_none=True)

    materials_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Material settings</h4>
    </div>"""
    )

    kpoints_distance_description = ipw.HTML(
        """<p>
        Similarly, the <b>protocol</b> also defines the k-points mesh density.
        Untick the box to override the default
    </p>"""
    )

    _DEFAULT_SPIN_TYPE = "NONE"
    _DEFAULT_KPOINTS_DISTANCE = 0.15
    _DEFAULT_ELECTRONIC_TYPE = "METAL"

    def __init__(self, **kwargs):

        self._spin_type = ipw.Dropdown(
            options=["NONE", "COLLINEAR"],
            value="NONE",
            description="Spin Type:",
            style={"description_width": "initial"},
        )
        self._spin_type.observe(self.set_spin_type_trait, "value")
        self.set_spin_type_trait()

        self._electronic_type = ipw.Dropdown(
            options=["METAL", "INSULATOR"],
            value="METAL",
            description="Electronic Type:",
            style={"description_width": "initial"},
        )
        self._electronic_type.observe(self.set_electronic_type_trait, "value")
        self.set_electronic_type_trait()

        # K-points distance.
        self._set_kpoints_distance_automatically = ipw.Checkbox(
            description="Use default k-points distance.",
            indent=False,
            value=True,
        )
        self._kpoints_distance = ipw.FloatText(
            value=self._DEFAULT_KPOINTS_DISTANCE,
            step=0.05,
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
                PseudoFamilySelector(),
                self.kpoints_distance_description,
                ipw.HBox(
                    [self._set_kpoints_distance_automatically, self._kpoints_distance]
                ),
                self.materials_title,
                self._spin_type,
                self._electronic_type
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs
        )

    def set_spin_type_trait(self, _=None):
        self.spin_type = SpinType[self._spin_type.value.upper()]

    def set_electronic_type_trait(self, _=None):
        self.electronic_type = ElectronicType[self._electronic_type.value]

    def set_kpoints_distance_trait(self, _=None):
        self.kpoints_distance = (
            self._DEFAULT_KPOINTS_DISTANCE
            if self._set_kpoints_distance_automatically.value
            else self._kpoints_distance.value
        )

class CodesConfig(ipw.VBox):

    codes_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Codes</h4></div>"""
    )

    codes_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        Here you can select the codes to use for running the calculations. The
        codes on the local machine are installed and selected by default, but
        you can set new ones for each of the codes by clicking on "Setup new
        code".</div>"""
    )

    def __init__(self, **kwargs):

        self.pw = CodeDropdown(
            input_plugin="quantumespresso.pw",
            description="pw.x:",
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
            description="dos.x",
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
            description="projwfc.x",
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
                self.codes_title,
                self.codes_help,
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
            ],
            layout=ipw.Layout(min_height="250px"),
        )

        self.tab.set_title(0, "Workchain")

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
            self.tab.set_title(1, "Advanced settings")
            self.tab.set_title(2, "Select codes")
            self.tab.set_title(3, "Compute resources")
            self.tab.children = [
                self.workchain_config,
                self.options_config,
                self.codes_selector,
                self.resources_config,
            ]
        else:
            self.tab.set_title(0, "Workchain")
            self.tab.children = [
                self.workchain_config,
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
        self.workchain_config.geo_opt_run.observe(update, ["value"])
        self.workchain_config.simulation_protocol.observe(update, ["value"])
        # "Extra" parameters
        self.workchain_config.bands_run.observe(update, ["value"])
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
                    if self.workchain_config.geo_opt_run.value
                    else RelaxType["NONE"],
                    spin_type=self.options_config.spin_type,
                    # "extra" parameters
                    run_bands=self.workchain_config.bands_run.value,
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
            self.workchain_config.geo_opt_run.value = (
                relax_type is not RelaxType["NONE"]
            )

            self.options_config.spin_type = bp.get("spin_type", SpinType["NONE"])
            # "extra" parameters
            self.workchain_config.bands_run.value = bp["run_bands"]
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
