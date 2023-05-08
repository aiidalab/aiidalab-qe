# -*- coding: utf-8 -*-
"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""
from __future__ import annotations

import os
import typing as t
from dataclasses import dataclass

import ipywidgets as ipw
import traitlets
from aiida.common import NotExistent
from aiida.engine import ProcessBuilderNamespace, ProcessState, submit
from aiida.orm import WorkChainNode, load_code, load_node
from aiida.plugins import DataFactory
from aiida_quantumespresso.common.types import ElectronicType, RelaxType, SpinType
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_widgets_base import (
    AiidaNodeViewWidget,
    ComputationalResourcesWidget,
    ProcessMonitor,
    ProcessNodesTreeWidget,
    WizardAppWidgetStep,
)
from IPython.display import display

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
from aiidalab_qe.app.pseudos import PseudoFamilySelector
from aiidalab_qe.app.setup_codes import QESetupWidget
from aiidalab_qe.app.sssp import SSSPInstallWidget
from aiidalab_qe.app.widgets import ParallelizationSettings, ResourceSelectionWidget
from aiidalab_qe.workflows import QeAppWorkChain

StructureData = DataFactory("core.structure")
Float = DataFactory("core.float")
Dict = DataFactory("core.dict")
Str = DataFactory("core.str")

PROTOCOL_PSEUDO_MAP = {
    "fast": "SSSP/1.2/PBE/efficiency",
    "moderate": "SSSP/1.2/PBE/efficiency",
    "precise": "SSSP/1.2/PBE/precision",
}


# The static input parameters for the QE App WorkChain
# The dataclass does not include codes and structure which will be set
# from widgets separately.
# Relax type, electronic type, spin type, are str because they are used also
# for serialized input of extras attributes of the workchain
@dataclass(frozen=True)
class QeWorkChainParameters:
    protocol: str
    relax_type: str
    properties: t.List[str]
    spin_type: str
    electronic_type: str
    overrides: t.Dict[str, t.Any]


class WorkChainSettings(ipw.VBox):
    structure_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Structure</h4></div>"""
    )
    structure_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
        You have three options:<br>
        (1) Structure as is: perform a self consistent calculation using the structure provided as input.<br>
        (2) Atomic positions: perform a full relaxation of the internal atomic coordinates. <br>
        (3) Full geometry: perform a full relaxation for both the internal atomic coordinates and the cell vectors. </div>"""
    )
    materials_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
        Below you can indicate both if the material should be treated as an insulator
        or a metal (if in doubt, choose "Metal"),
        and if it should be studied with magnetization/spin polarization,
        switch magnetism On or Off (On is at least twice more costly).
        </div>"""
    )

    properties_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Properties</h4></div>"""
    )
    properties_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px">
        The band structure workflow will
        automatically detect the default path in reciprocal space using the
        <a href="https://www.materialscloud.org/work/tools/seekpath" target="_blank">
        SeeK-path tool</a>.</div>"""
    )

    protocol_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Protocol</h4></div>"""
    )
    protocol_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 6px; padding-bottom: 0px">
        The "moderate" protocol represents a trade-off between
        accuracy and speed. Choose the "fast" protocol for a faster calculation
        with less precision and the "precise" protocol to aim at best accuracy (at the price of longer/costlier calculations).</div>"""
    )

    def __init__(self, **kwargs):
        # RelaxType: degrees of freedom in geometry optimization
        self.relax_type = ipw.ToggleButtons(
            options=[
                ("Structure as is", "none"),
                ("Atomic positions", "positions"),
                ("Full geometry", "positions_cell"),
            ],
            value="positions_cell",
        )

        # SpinType: magnetic properties of material
        self.spin_type = ipw.ToggleButtons(
            options=[("Off", "none"), ("On", "collinear")],
            value=DEFAULT_PARAMETERS["spin_type"],
            style={"description_width": "initial"},
        )

        # ElectronicType: electronic properties of material
        self.electronic_type = ipw.ToggleButtons(
            options=[("Metal", "metal"), ("Insulator", "insulator")],
            value=DEFAULT_PARAMETERS["electronic_type"],
            style={"description_width": "initial"},
        )

        # Checkbox to see if the band structure should be calculated
        self.bands_run = ipw.Checkbox(
            description="",
            indent=False,
            value=True,
            layout=ipw.Layout(max_width="10%"),
        )

        # Checkbox to see if the PDOS should be calculated
        self.pdos_run = ipw.Checkbox(
            description="",
            indent=False,
            value=True,
            layout=ipw.Layout(max_width="10%"),
        )

        # Work chain protocol
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )
        super().__init__(
            children=[
                self.structure_title,
                self.structure_help,
                self.relax_type,
                self.materials_help,
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "Electronic Type:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="120px"
                            ),
                        ),
                        self.electronic_type,
                    ]
                ),
                ipw.HBox(
                    children=[
                        ipw.Label(
                            "Magnetism:",
                            layout=ipw.Layout(
                                justify_content="flex-start", width="120px"
                            ),
                        ),
                        self.spin_type,
                    ]
                ),
                self.properties_title,
                ipw.HTML("Select which properties to calculate:"),
                ipw.HBox(children=[ipw.HTML("<b>Band structure</b>"), self.bands_run]),
                ipw.HBox(
                    children=[
                        ipw.HTML("<b>Projected density of states</b>"),
                        self.pdos_run,
                    ]
                ),
                self.properties_help,
                self.protocol_title,
                ipw.HTML("Select the protocol:", layout=ipw.Layout(flex="1 1 auto")),
                self.workchain_protocol,
                self.protocol_help,
            ],
            **kwargs,
        )

    def _update_settings(self, **kwargs):
        """Update the settings based on the given dict."""
        for key in [
            "relax_type",
            "spin_type",
            "electronic_type",
            "bands_run",
            "pdos_run",
            "workchain_protocol",
        ]:
            if key in kwargs:
                getattr(self, key).value = kwargs[key]


class SmearingSettings(ipw.VBox):
    smearing_description = ipw.HTML(
        """<p>
        The smearing type and width is set by the chosen <b>protocol</b>.
        Tick the box to override the default, not advised unless you've mastered <b>smearing effects</b> (click <a href="http://theossrv1.epfl.ch/Main/ElectronicTemperature"
        target="_blank">here</a> for a discussion).
    </p>"""
    )

    # The default of `smearing` and `degauss` the type and width
    # must be linked to the `protocol`
    degauss_default = traitlets.Float(default_value=0.01)
    smearing_default = traitlets.Unicode(default_value="cold")

    def __init__(self, **kwargs):
        self.override_protocol_smearing = ipw.Checkbox(
            description="Override",
            indent=False,
            value=False,
        )
        self.smearing = ipw.Dropdown(
            options=["cold", "gaussian", "fermi-dirac", "methfessel-paxton"],
            value=self.smearing_default,
            description="Smearing type:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.degauss = ipw.FloatText(
            value=self.degauss_default,
            step=0.005,
            description="Smearing width (Ry):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override_protocol_smearing, "value"),
            (self.degauss, "disabled"),
            lambda override: not override,
        )
        ipw.dlink(
            (self.override_protocol_smearing, "value"),
            (self.smearing, "disabled"),
            lambda override: not override,
        )
        self.degauss.observe(self.set_smearing, "value")
        self.smearing.observe(self.set_smearing, "value")
        self.override_protocol_smearing.observe(self.set_smearing, "value")

        super().__init__(
            children=[
                self.smearing_description,
                ipw.HBox(
                    [self.override_protocol_smearing, self.smearing, self.degauss]
                ),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

    def set_smearing(self, _=None):
        self.degauss.value = (
            self.degauss.value
            if self.override_protocol_smearing.value
            else self.degauss_default
        )
        self.smearing.value = (
            self.smearing.value
            if self.override_protocol_smearing.value
            else self.smearing_default
        )

    def _update_settings(self, **kwargs):
        """Update the smearing and degauss values by the given keyword arguments
        This is the same as the `set_smearing` method but without the observer.
        Therefore the override checkbox is not updated and defaults to True"""
        self.override_protocol_smearing.value = True

        with self.hold_trait_notifications():
            if "smearing" in kwargs:
                self.smearing.value = kwargs["smearing"]

            if "degauss" in kwargs:
                self.degauss.value = kwargs["degauss"]


class KpointSettings(ipw.VBox):
    kpoints_distance_description = ipw.HTML(
        """<div>
        The k-points mesh density of the SCF calculation is set by the <b>protocol</b>.
        The value below represents the maximum distance between the k-points in each direction of reciprocal space.
        Tick the box to override the default, smaller is more accurate and costly. </div>"""
    )

    # The default of `kpoints_distance` must be linked to the `protocol`
    kpoints_distance_default = traitlets.Float(default_value=0.15)

    def __init__(self, **kwargs):
        self.override_protocol_kpoints = ipw.Checkbox(
            description="Override",
            indent=False,
            value=False,
        )
        self.kpoints_distance = ipw.FloatText(
            value=self.kpoints_distance_default,
            step=0.05,
            description="K-points distance (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.dlink(
            (self.override_protocol_kpoints, "value"),
            (self.kpoints_distance, "disabled"),
            lambda override: not override,
        )
        self.kpoints_distance.observe(self.set_kpoints_distance, "value")
        self.override_protocol_kpoints.observe(self.set_kpoints_distance, "value")
        self.observe(self.set_kpoints_distance, "kpoints_distance_default")

        super().__init__(
            children=[
                self.kpoints_distance_description,
                ipw.HBox([self.override_protocol_kpoints, self.kpoints_distance]),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )

    def set_kpoints_distance(self, _=None):
        self.kpoints_distance.value = (
            self.kpoints_distance.value
            if self.override_protocol_kpoints.value
            else self.kpoints_distance_default
        )

    def _update_settings(self, **kwargs):
        """Update the kpoints_distance value by the given keyword arguments.
        This is the same as the `set_kpoints_distance` method but without the observer.
        """
        self.override_protocol_kpoints.value = True
        if "kpoints_distance" in kwargs:
            self.kpoints_distance.value = kwargs["kpoints_distance"]


class ConfigureQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    confirmed = traitlets.Bool()
    previous_step_state = traitlets.UseEnum(WizardAppWidgetStep.State)
    workchain_settings = traitlets.Instance(WorkChainSettings, allow_none=True)
    kpoints_settings = traitlets.Instance(KpointSettings, allow_none=True)
    smearing_settings = traitlets.Instance(SmearingSettings, allow_none=True)
    pseudo_family_selector = traitlets.Instance(PseudoFamilySelector, allow_none=True)

    def __init__(self, **kwargs):
        self.workchain_settings = WorkChainSettings()
        self.workchain_settings.relax_type.observe(self._update_state, "value")
        self.workchain_settings.bands_run.observe(self._update_state, "value")
        self.workchain_settings.pdos_run.observe(self._update_state, "value")

        self.kpoints_settings = KpointSettings()
        self.smearing_settings = SmearingSettings()
        self.pseudo_family_selector = PseudoFamilySelector()

        ipw.dlink(
            (self.workchain_settings.workchain_protocol, "value"),
            (self.kpoints_settings, "kpoints_distance_default"),
            lambda protocol: PwBaseWorkChain.get_protocol_inputs(protocol)[
                "kpoints_distance"
            ],
        )

        ipw.dlink(
            (self.workchain_settings.workchain_protocol, "value"),
            (self.smearing_settings, "degauss_default"),
            lambda protocol: PwBaseWorkChain.get_protocol_inputs(protocol)["pw"][
                "parameters"
            ]["SYSTEM"]["degauss"],
        )

        ipw.dlink(
            (self.workchain_settings.workchain_protocol, "value"),
            (self.smearing_settings, "smearing_default"),
            lambda protocol: PwBaseWorkChain.get_protocol_inputs(protocol)["pw"][
                "parameters"
            ]["SYSTEM"]["smearing"],
        )

        self.tab = ipw.Tab(
            children=[
                self.workchain_settings,
                ipw.VBox(
                    children=[
                        self.pseudo_family_selector,
                        self.kpoints_settings,
                        self.smearing_settings,
                    ]
                ),
            ],
            layout=ipw.Layout(min_height="250px"),
        )

        self.tab.set_title(0, "Workflow")
        self.tab.set_title(1, "Advanced settings")

        self._submission_blocker_messages = ipw.HTML()

        self.confirm_button = ipw.Button(
            description="Confirm",
            tooltip="Confirm the currently selected settings and go to the next step.",
            button_style="success",
            icon="check-circle",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )

        self.confirm_button.on_click(self.confirm)

        super().__init__(
            children=[
                self.tab,
                self._submission_blocker_messages,
                self.confirm_button,
            ],
            **kwargs,
        )

    @traitlets.observe("previous_step_state")
    def _observe_previous_step_state(self, change):
        self._update_state()

    def set_input_parameters(self, parameters):
        """Set the inputs in the GUI based on a set of parameters."""

        with self.hold_trait_notifications():
            # Work chain settings
            self.workchain_settings.relax_type.value = parameters["relax_type"]
            self.workchain_settings.spin_type.value = parameters["spin_type"]
            self.workchain_settings.electronic_type.value = parameters[
                "electronic_type"
            ]
            self.workchain_settings.bands_run.value = parameters["run_bands"]
            self.workchain_settings.pdos_run.value = parameters["run_pdos"]
            self.workchain_settings.workchain_protocol.value = parameters["protocol"]

            # Advanced settings
            self.pseudo_family_selector.value = parameters["pseudo_family"]
            if parameters.get("kpoints_distance_override", None) is not None:
                self.kpoints_settings.kpoints_distance.value = parameters[
                    "kpoints_distance_override"
                ]
                self.kpoints_settings.override_protocol_kpoints.value = True
            if parameters.get("degauss_override", None) is not None:
                self.smearing_settings.degauss.value = parameters["degauss_override"]
                self.smearing_settings.override_protocol_smearing.value = True
            if parameters.get("smearing_override", None) is not None:
                self.smearing_settings.smearing.value = parameters["smearing_override"]
                self.smearing_settings.override_protocol_smearing.value = True

    def _update_state(self, _=None):
        if self.previous_step_state == self.State.SUCCESS:
            self.confirm_button.disabled = False
            self._submission_blocker_messages.value = ""
            self.state = self.State.CONFIGURED
        elif self.previous_step_state == self.State.FAIL:
            self.state = self.State.FAIL
        else:
            self.confirm_button.disabled = True
            self.state = self.State.INIT
            self.set_input_parameters(DEFAULT_PARAMETERS)

    def confirm(self, _=None):
        self.confirm_button.disabled = False
        self.state = self.State.SUCCESS

    @traitlets.default("state")
    def _default_state(self):
        return self.State.INIT

    def reset(self):
        with self.hold_trait_notifications():
            self.set_input_parameters(DEFAULT_PARAMETERS)


class SubmitQeAppWorkChainStep(ipw.VBox, WizardAppWidgetStep):
    """Step for submission of a bands workchain."""

    codes_title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 0px">
        <h4>Codes</h4></div>"""
    )
    codes_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom:
        10px"> Select the code to use for running the calculations. The codes
        on the local machine (localhost) are installed by default, but you can
        configure new ones on potentially more powerful machines by clicking on
        "Setup new code".</div>"""
    )

    # This number provides a rough estimate for how many MPI tasks are needed
    # for a given structure.
    NUM_SITES_PER_MPI_TASK_DEFAULT = 6

    # Warn the user if they are trying to run calculations for a large
    # structure on localhost.
    RUN_ON_LOCALHOST_NUM_SITES_WARN_THRESHOLD = 10

    # Put a limit on how many MPI tasks you want to run per k-pool by default
    MAX_MPI_PER_POOL = 20

    input_structure = traitlets.Instance(StructureData, allow_none=True)
    process = traitlets.Instance(WorkChainNode, allow_none=True)
    previous_step_state = traitlets.UseEnum(WizardAppWidgetStep.State)
    workchain_settings = traitlets.Instance(WorkChainSettings, allow_none=True)
    kpoints_settings = traitlets.Instance(KpointSettings, allow_none=True)
    smearing_settings = traitlets.Instance(SmearingSettings, allow_none=True)
    pseudo_family_selector = traitlets.Instance(PseudoFamilySelector, allow_none=True)
    _submission_blockers = traitlets.List(traitlets.Unicode())

    def __init__(self, qe_auto_setup=True, **kwargs):
        self.message_area = ipw.Output()
        self._submission_blocker_messages = ipw.HTML()

        self.pw_code = ComputationalResourcesWidget(
            description="pw.x:", default_calc_job_plugin="quantumespresso.pw"
        )
        self.dos_code = ComputationalResourcesWidget(
            description="dos.x:",
            default_calc_job_plugin="quantumespresso.dos",
        )
        self.projwfc_code = ComputationalResourcesWidget(
            description="projwfc.x:",
            default_calc_job_plugin="quantumespresso.projwfc",
        )

        self.resources_config = ResourceSelectionWidget()
        self.parallelization = ParallelizationSettings()

        self.set_selected_codes(DEFAULT_PARAMETERS)
        self.set_resource_defaults()

        self.pw_code.observe(self._update_state, "value")
        self.pw_code.observe(self._update_resources, "value")
        self.dos_code.observe(self._update_state, "value")
        self.projwfc_code.observe(self._update_state, "value")

        self.submit_button = ipw.Button(
            description="Submit",
            tooltip="Submit the calculation with the selected parameters.",
            icon="play",
            button_style="success",
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
            disabled=True,
        )

        self.submit_button.on_click(self._on_submit_button_clicked)

        # The SSSP installation status widget shows the installation status of
        # the SSSP pseudo potentials and triggers the installation in case that
        # they are not yet installed. The widget will remain in a "busy" state
        # in case that the installation was already triggered elsewhere, e.g.,
        # by the start up scripts.  The submission is blocked while the
        # potentials are not yet installed.
        self.sssp_installation_status = SSSPInstallWidget(auto_start=qe_auto_setup)
        self.sssp_installation_status.observe(self._update_state, ["busy", "installed"])
        self.sssp_installation_status.observe(self._toggle_install_widgets, "installed")

        # The QE setup widget checks whether there are codes that match specific
        # expected labels (e.g. "pw-7.0@localhost") and triggers both the
        # installation of QE into a dedicated conda environment and the setup of
        # the codes in case that they are not already configured.
        self.qe_setup_status = QESetupWidget(auto_start=qe_auto_setup)
        self.qe_setup_status.observe(self._update_state, "busy")
        self.qe_setup_status.observe(self._toggle_install_widgets, "installed")
        self.qe_setup_status.observe(self._auto_select_code, "installed")

        super().__init__(
            children=[
                self.codes_title,
                self.codes_help,
                self.pw_code,
                self.dos_code,
                self.projwfc_code,
                self.resources_config,
                self.parallelization,
                self.message_area,
                self.sssp_installation_status,
                self.qe_setup_status,
                self._submission_blocker_messages,
                self.submit_button,
            ]
        )

    @traitlets.observe("_submission_blockers")
    def _observe_submission_blockers(self, change):
        if change["new"]:
            fmt_list = "\n".join((f"<li>{item}</li>" for item in sorted(change["new"])))
            self._submission_blocker_messages.value = f"""
                <div class="alert alert-info">
                <strong>The submission is blocked, due to the following reason(s):</strong>
                <ul>{fmt_list}</ul></div>"""
        else:
            self._submission_blocker_messages.value = ""

    def _identify_submission_blockers(self):
        # Do not submit while any of the background setup processes are running.
        if self.qe_setup_status.busy or self.sssp_installation_status.busy:
            yield "Background setup processes must finish."

        # No code selected (this is ignored while the setup process is running).
        if self.pw_code.value is None and not self.qe_setup_status.busy:
            yield ("No pw code selected")

        # No code selected for pdos (this is ignored while the setup process is running).
        if (
            self.workchain_settings.pdos_run.value
            and (self.dos_code.value is None or self.projwfc_code.value is None)
            and not self.qe_setup_status.busy
        ):
            yield "Calculating the PDOS requires both dos.x and projwfc.x to be set."

        # SSSP library not installed
        if not self.sssp_installation_status.installed:
            yield "The SSSP library is not installed."

        if (
            self.workchain_settings.pdos_run.value
            and not any(
                [
                    self.pw_code.value is None,
                    self.dos_code.value is None,
                    self.projwfc_code.value is None,
                ]
            )
            and len(
                set(
                    (
                        load_code(self.pw_code.value).computer.pk,
                        load_code(self.dos_code.value).computer.pk,
                        load_code(self.projwfc_code.value).computer.pk,
                    )
                )
            )
            != 1
        ):
            yield (
                "All selected codes must be installed on the same computer. This is because the "
                "PDOS calculations rely on large files that are not retrieved by AiiDA."
            )

    def _update_state(self, _=None):
        # If the previous step has failed, this should fail as well.
        if self.previous_step_state is self.State.FAIL:
            self.state = self.State.FAIL
            return
        # Do not interact with the user if they haven't successfully completed the previous step.
        elif self.previous_step_state is not self.State.SUCCESS:
            self.state = self.State.INIT
            return

        # Process is already running.
        if self.process is not None:
            self.state = self.State.SUCCESS
            return

        blockers = list(self._identify_submission_blockers())
        if any(blockers):
            self._submission_blockers = blockers
            self.state = self.State.READY
            return

        self._submission_blockers = []
        self.state = self.state.CONFIGURED

    def _toggle_install_widgets(self, change):
        if change["new"]:
            self.children = [
                child for child in self.children if child is not change["owner"]
            ]

    def _auto_select_code(self, change):
        if change["new"] and not change["old"]:
            for code in [
                "pw_code",
                "dos_code",
                "projwfc_code",
            ]:
                try:
                    code_widget = getattr(self, code)
                    code_widget.refresh()
                    code_widget.value = load_code(DEFAULT_PARAMETERS[code]).uuid
                except NotExistent:
                    pass

    _ALERT_MESSAGE = """
        <div class="alert alert-{alert_class} alert-dismissible">
        <a href="#" class="close" data-dismiss="alert" aria-label="close">&times;</a>
        <span class="closebtn" onclick="this.parentElement.style.display='none';">&times;</span>
        <strong>{message}</strong>
        </div>"""

    def _show_alert_message(self, message, alert_class="info"):
        with self.message_area:
            display(
                ipw.HTML(
                    self._ALERT_MESSAGE.format(alert_class=alert_class, message=message)
                )
            )

    def _update_resources(self, change):
        if change["new"] and (
            change["old"] is None
            or load_code(change["new"]).computer.pk
            != load_code(change["old"]).computer.pk
        ):
            self.set_resource_defaults(load_code(change["new"]).computer)

    def set_resource_defaults(self, computer=None):
        if computer is None or computer.hostname == "localhost":
            self.resources_config.num_nodes.disabled = True
            self.resources_config.num_nodes.value = 1
            self.resources_config.num_cpus.max = os.cpu_count()
            self.resources_config.num_cpus.value = 1
            self.resources_config.num_cpus.description = "CPUs"
            self.parallelization.npools.value = 1
        else:
            default_mpiprocs = computer.get_default_mpiprocs_per_machine()
            self.resources_config.num_nodes.disabled = False
            self.resources_config.num_cpus.max = default_mpiprocs
            self.resources_config.num_cpus.value = default_mpiprocs
            self.resources_config.num_cpus.description = "CPUs/node"
            self.parallelization.npools.value = self._get_default_parallelization()

        self._check_resources()

    def _get_default_parallelization(self):
        """A _very_ rudimentary approach for obtaining a minimal npools setting."""
        num_mpiprocs = (
            self.resources_config.num_nodes.value * self.resources_config.num_cpus.value
        )

        for i in range(1, num_mpiprocs + 1):
            if num_mpiprocs % i == 0 and num_mpiprocs // i < self.MAX_MPI_PER_POOL:
                return i

    def _check_resources(self):
        """Check whether the currently selected resources will be sufficient and warn if not."""
        if not self.pw_code.value:
            return  # No code selected, nothing to do.

        num_cpus = self.resources_config.num_cpus.value
        on_localhost = load_node(self.pw_code.value).computer.hostname == "localhost"
        if self.pw_code.value and on_localhost and num_cpus > 1:
            self._show_alert_message(
                "The selected code would be executed on the local host, but "
                "the number of CPUs is larger than one. Please review "
                "the configuration and consider to select a code that runs "
                "on a larger system if necessary.",
                alert_class="warning",
            )
        elif (
            self.input_structure
            and on_localhost
            and len(self.input_structure.sites)
            > self.RUN_ON_LOCALHOST_NUM_SITES_WARN_THRESHOLD
        ):
            self._show_alert_message(
                "The selected code would be executed on the local host, but the "
                "number of sites of the selected structure is relatively large. "
                "Consider to select a code that runs on a larger system if "
                "necessary.",
                alert_class="warning",
            )

    @traitlets.observe("state")
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            self.submit_button.disabled = change["new"] != self.State.CONFIGURED

    @traitlets.observe("previous_step_state")
    def _observe_input_structure(self, _):
        self._update_state()
        self.set_pdos_status()

    @traitlets.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            process_node = change["new"]
            if process_node is not None:
                self.input_structure = process_node.inputs.structure
                builder_parameters = process_node.base.extras.get(
                    "builder_parameters", None
                )
                if builder_parameters is not None:
                    self.set_selected_codes(builder_parameters)
            self._update_state()

    def _on_submit_button_clicked(self, _):
        self.submit_button.disabled = True
        self.submit()

    def set_selected_codes(self, parameters):
        """Set the inputs in the GUI based on a set of parameters."""

        # Codes
        def _get_code_uuid(code):
            if code is not None:
                try:
                    return load_code(code).uuid
                except NotExistent:
                    return None

        with self.hold_trait_notifications():
            # Codes
            self.pw_code.value = _get_code_uuid(parameters["pw_code"])
            self.dos_code.value = _get_code_uuid(parameters["dos_code"])
            self.projwfc_code.value = _get_code_uuid(parameters["projwfc_code"])

    def set_pdos_status(self):
        if self.workchain_settings.pdos_run.value:
            self.dos_code.code_select_dropdown.disabled = False
            self.projwfc_code.code_select_dropdown.disabled = False
        else:
            self.dos_code.code_select_dropdown.disabled = True
            self.projwfc_code.code_select_dropdown.disabled = True

    def submit(self, _=None):
        """Submit the work chain with the current inputs."""
        builder = self._create_builder()
        extra_parameters = self._create_extra_report_parameters()

        with self.hold_trait_notifications():
            self.process = submit(builder)

            # Set the builder parameters on the work chain
            builder_parameters = self._extract_report_parameters(
                builder, extra_parameters
            )
            self.process.base.extras.set("builder_parameters", builder_parameters)

    def _get_qe_workchain_parameters(self) -> QeWorkChainParameters:
        """Get the parameters of the `QeWorkChain` from widgets."""
        # create the override parameters for sub PwBaseWorkChain
        pw_overrides = {"base": {}, "scf": {}, "nscf": {}, "band": {}}
        for key in ["base", "scf", "nscf", "band"]:
            if self.pseudo_family_selector.override_protocol_pseudo_family.value:
                pw_overrides[key]["pseudo_family"] = self.pseudo_family_selector.value

        for key in ["base", "scf"]:
            if self.kpoints_settings.override_protocol_kpoints.value:
                pw_overrides[key][
                    "kpoints_distance"
                ] = self.kpoints_settings.kpoints_distance.value

            if (
                self.smearing_settings.override_protocol_smearing.value
                and self.workchain_settings.electronic_type.value == "metal"
            ):
                # smearing type setting
                pw_overrides[key].setdefault("pw", {}).setdefault(
                    "parameters", {}
                ).setdefault("SYSTEM", {})[
                    "smearing"
                ] = self.smearing_settings.smearing.value

                # smearing degauss setting
                pw_overrides[key].setdefault("pw", {}).setdefault(
                    "parameters", {}
                ).setdefault("SYSTEM", {})[
                    "degauss"
                ] = self.smearing_settings.degauss.value

        overrides = {
            "relax": {
                "base": pw_overrides["base"],
            },
            "bands": {
                "scf": pw_overrides["scf"],
                "bands": pw_overrides["band"],
            },
            "pdos": {
                "scf": pw_overrides["scf"],
                "nscf": pw_overrides["nscf"],
            },
        }

        # Work chain settings
        relax_type = self.workchain_settings.relax_type.value
        electronic_type = self.workchain_settings.electronic_type.value
        spin_type = self.workchain_settings.spin_type.value

        run_bands = self.workchain_settings.bands_run.value
        run_pdos = self.workchain_settings.pdos_run.value
        protocol = self.workchain_settings.workchain_protocol.value

        properties = []

        if run_bands:
            properties.append("bands")
        if run_pdos:
            properties.append("pdos")

        if RelaxType(relax_type) is not RelaxType.NONE or not (run_bands or run_pdos):
            properties.append("relax")

        return QeWorkChainParameters(
            protocol=protocol,
            relax_type=relax_type,
            properties=properties,
            spin_type=spin_type,
            electronic_type=electronic_type,
            overrides=overrides,
        )

    def _create_builder(self) -> ProcessBuilderNamespace:
        """Create the builder for the `QeAppWorkChain` submit."""
        pw_code = self.pw_code.value
        dos_code = self.dos_code.value
        projwfc_code = self.projwfc_code.value

        parameters = self._get_qe_workchain_parameters()

        builder = QeAppWorkChain.get_builder_from_protocol(
            structure=self.input_structure,
            pw_code=load_code(pw_code),
            dos_code=load_code(dos_code),
            projwfc_code=load_code(projwfc_code),
            protocol=parameters.protocol,
            relax_type=RelaxType(parameters.relax_type),
            properties=parameters.properties,
            spin_type=SpinType(parameters.spin_type),
            electronic_type=ElectronicType(parameters.electronic_type),
            overrides=parameters.overrides,
        )

        resources = {
            "num_machines": self.resources_config.num_nodes.value,
            "num_mpiprocs_per_machine": self.resources_config.num_cpus.value,
        }

        npool = self.parallelization.npools.value
        self._update_builder(builder, resources, npool, self.MAX_MPI_PER_POOL)

        return builder

    def _update_builder(self, buildy, resources, npools, max_mpi_per_pool):
        """Update the resources and parallelization of the ``QeAppWorkChain`` builder."""
        for k, v in buildy.items():
            if isinstance(v, (dict, ProcessBuilderNamespace)):
                if k == "pw" and v["pseudos"]:
                    v["parallelization"] = Dict(dict={"npool": npools})
                if k == "projwfc":
                    v["settings"] = Dict(dict={"cmdline": ["-nk", str(npools)]})
                if k == "dos":
                    v["metadata"]["options"]["resources"] = {
                        "num_machines": 1,
                        "num_mpiprocs_per_machine": min(
                            max_mpi_per_pool,
                            resources["num_mpiprocs_per_machine"],
                        ),
                    }
                    # Continue to the next item to avoid overriding the resources in the
                    # recursive `update_builder` call.
                    continue
                if k == "resources":
                    buildy["resources"] = resources
                else:
                    self._update_builder(v, resources, npools, max_mpi_per_pool)

    def _create_extra_report_parameters(self) -> dict[str, t.Any]:
        """This method will also create a dictionary of the parameters that were not
        readably represented in the builder, which will be used to the report.
        It is stored in the `extra_report_parameters`.
        """
        qe_workchain_parameters = self._get_qe_workchain_parameters()

        # Construct the extra report parameters needed for the report
        extra_report_parameters = {
            "relax_type": qe_workchain_parameters.relax_type,
            "electronic_type": qe_workchain_parameters.electronic_type,
            "spin_type": qe_workchain_parameters.spin_type,
            "protocol": qe_workchain_parameters.protocol,
        }

        # update pseudo family information to extra_report_parameters
        if self.pseudo_family_selector.override_protocol_pseudo_family.value:
            # If the pseudo family is overridden, use that
            pseudo_family = self.pseudo_family_selector.value
        else:
            # otherwise extract the information from protocol
            pseudo_family = PROTOCOL_PSEUDO_MAP[qe_workchain_parameters.protocol]

        pseudo_family_info = pseudo_family.split("/")
        extra_report_parameters.update(
            {
                "pseudo_family": pseudo_family,
                "pseudo_library": pseudo_family_info[0],
                "pseudo_version": pseudo_family_info[1],
                "functional": pseudo_family_info[2],
                "pseudo_protocol": pseudo_family_info[3],
            }
        )

        # store codes info into extra_report_parameters for loading the process
        pw_code = self.pw_code.value
        dos_code = self.dos_code.value
        projwfc_code = self.projwfc_code.value

        extra_report_parameters.update(
            {
                "pw_code": pw_code,
                "dos_code": dos_code,
                "projwfc_code": projwfc_code,
            }
        )

        return extra_report_parameters

    @staticmethod
    def _extract_report_parameters(
        builder, extra_report_parameters
    ) -> dict[str, t.Any]:
        """Extract (recover) the parameters for report from the builder.

        There are some parameters that are not stored in the builder, but can be extracted
        directly from the widgets, such as the ``pseudo_family`` and ``relax_type``.
        """
        parameters = {
            "run_relax": "relax" in builder.properties,
            "run_bands": "bands" in builder.properties,
            "run_pdos": "pdos" in builder.properties,
        }

        # Extract the pw calculation parameters from the builder

        # energy_cutoff is same for all pw calculations when pseudopotentials are fixed
        # as well as the smearing settings (semaring and degauss) and scf kpoints distance
        # read from the first pw calculation of relax workflow.
        # It is safe then to extract these parameters from the first pw calculation, since the
        # builder is anyway set with subworkchain inputs even it is not run which controlled by
        # the properties inputs.
        energy_cutoff_wfc = builder.relax.base["pw"]["parameters"]["SYSTEM"]["ecutwfc"]
        energy_cutoff_rho = builder.relax.base["pw"]["parameters"]["SYSTEM"]["ecutrho"]
        occupation = builder.relax.base["pw"]["parameters"]["SYSTEM"]["occupations"]
        scf_kpoints_distance = builder.relax.base.kpoints_distance.value

        parameters.update(
            {
                "energy_cutoff_wfc": energy_cutoff_wfc,
                "energy_cutoff_rho": energy_cutoff_rho,
                "occupation": occupation,
                "scf_kpoints_distance": scf_kpoints_distance,
            }
        )

        if occupation == "smearing":
            parameters["degauss"] = builder.relax.base["pw"]["parameters"]["SYSTEM"][
                "degauss"
            ]
            parameters["smearing"] = builder.relax.base["pw"]["parameters"]["SYSTEM"][
                "smearing"
            ]

        parameters[
            "bands_kpoints_distance"
        ] = builder.bands.bands_kpoints_distance.value
        parameters["nscf_kpoints_distance"] = builder.pdos.nscf.kpoints_distance.value

        # parameters from extra_report_parameters
        for k, v in extra_report_parameters.items():
            parameters.update({k: v})

        return parameters

    def reset(self):
        with self.hold_trait_notifications():
            self.process = None
            self.input_structure = None


class ViewQeAppWorkChainStatusAndResultsStep(ipw.VBox, WizardAppWidgetStep):
    process = traitlets.Unicode(allow_none=True)

    def __init__(self, **kwargs):
        self.process_tree = ProcessNodesTreeWidget()
        ipw.dlink(
            (self, "process"),
            (self.process_tree, "value"),
        )

        self.node_view = AiidaNodeViewWidget(layout={"width": "auto", "height": "auto"})
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
        ipw.dlink((self, "process"), (self.process_monitor, "value"))

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
            process = load_node(self.process)
            process_state = process.process_state
            if process_state in (
                ProcessState.CREATED,
                ProcessState.RUNNING,
                ProcessState.WAITING,
            ):
                self.state = self.State.ACTIVE
            elif (
                process_state in (ProcessState.EXCEPTED, ProcessState.KILLED)
                or process.is_failed
            ):
                self.state = self.State.FAIL
            elif process.is_finished_ok:
                self.state = self.State.SUCCESS

    @traitlets.observe("process")
    def _observe_process(self, change):
        self._update_state()
