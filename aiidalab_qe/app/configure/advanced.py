import ipywidgets as ipw
from aiida.plugins import DataFactory
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from IPython.display import clear_output, display

from aiidalab_qe.app.common.panel import Panel
from aiidalab_qe.app.configure.pseudos import PseudoFamilySelector
from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS as QEAPP_DEFAULT_PARAMETERS

StructureData = DataFactory("core.structure")

DEFAULT_PARAMETERS = QEAPP_DEFAULT_PARAMETERS["advanced"]


class AdvancedSettings(Panel):
    title = "Advanced Settings"

    subtitle = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Advanced Settings</h4></div>"""
    )
    subdescription = ipw.HTML(
        """Select the advanced settings for the <b>pw.x</b> code. Tick the box to override the default."""
    )

    smearing_description = ipw.HTML(
        """<p>
        The smearing type and width is set by the chosen <b>protocol</b>.
        About <b>smearing effects</b> (click <a href="http://theossrv1.epfl.ch/Main/ElectronicTemperature"
        target="_blank">here</a> for a discussion).
    </p>"""
    )
    kpoints_distance_description = ipw.HTML(
        """<div>
        The k-points mesh density of the SCF calculation is set by the <b>protocol</b>.
        The value below represents the maximum distance between the k-points in each direction of reciprocal space.
        Smaller is more accurate and costly. </div>"""
    )
    tot_charge_description = ipw.HTML("""<p></p>""")

    def __init__(self, **kwargs):
        self.override = ipw.Checkbox(
            description="Override",
            indent=False,
            value=False,
        )
        # Work chain protocol
        self.workchain_protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
            value="moderate",
        )
        self.pseudo_family_selector = PseudoFamilySelector()
        #
        self.smearing = ipw.Dropdown(
            options=["cold", "gaussian", "fermi-dirac", "methfessel-paxton"],
            value="cold",
            description="<b>Smearing type</b>:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.degauss = ipw.FloatText(
            value=0.01,
            step=0.005,
            description="<b>Smearing width</b> (Ry):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.kpoints_distance = ipw.BoundedFloatText(
            value=0.15,
            min=0,
            step=0.05,
            description="<b>K-points distance</b> (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.tot_charge = ipw.BoundedFloatText(
            value=QEAPP_DEFAULT_PARAMETERS["advanced"]["tot_charge"],
            min=-3,
            max=3,
            step=0.01,
            disabled=False,
            description="<b>Total charge</b>:",
            style={"description_width": "initial"},
        )
        self.magnetization = MagnetizationSettings()
        # override
        for item in [
            self.tot_charge,
            self.kpoints_distance,
            self.smearing,
            self.degauss,
        ]:
            ipw.dlink(
                (self.override, "value"),
                (item, "disabled"),
                lambda override: not override,
            )
        ipw.dlink(
            (self.override, "value"),
            (self.magnetization.override, "value"),
        )
        # update settings based on protocol
        ipw.dlink(
            (self.workchain_protocol, "value"),
            (self.kpoints_distance, "value"),
            lambda protocol: PwBaseWorkChain.get_protocol_inputs(protocol)[
                "kpoints_distance"
            ],
        )
        ipw.dlink(
            (self.workchain_protocol, "value"),
            (self.degauss, "value"),
            lambda protocol: PwBaseWorkChain.get_protocol_inputs(protocol)["pw"][
                "parameters"
            ]["SYSTEM"]["degauss"],
        )

        ipw.dlink(
            (self.workchain_protocol, "value"),
            (self.smearing, "value"),
            lambda protocol: PwBaseWorkChain.get_protocol_inputs(protocol)["pw"][
                "parameters"
            ]["SYSTEM"]["smearing"],
        )
        #
        self.children = [
            self.subtitle,
            ipw.HBox(
                [
                    self.subdescription,
                    self.override,
                ],
                layout=ipw.Layout(height="50px", justify_content="space-between"),
            ),
            # self.tot_charge_description,
            self.tot_charge,
            self.magnetization,
            self.kpoints_distance_description,
            self.kpoints_distance,
            self.smearing_description,
            self.smearing,
            self.degauss,
            self.pseudo_family_selector,
        ]
        super().__init__(
            **kwargs,
        )

    def get_panel_value(self):
        """Return the value of all the widgets in the panel as a dictionary.

        :return: a dictionary of the values of all the widgets in the panel.
        """
        # TODO insulator, no smearing, the occupation type is set to fixed, smearing and degauss should not be set
        # TODO "fixed occupations and lsda need tot_magnetization"
        if (
            self.parent is not None
            and self.parent.basic_settings.spin_type.value == "collinear"
        ):
            initial_magnetic_moments = self.magnetization.get_magnetization()
        else:
            initial_magnetic_moments = None
        parameters = {
            "pseudo_family": self.pseudo_family_selector.value,
            "kpoints_distance": self.kpoints_distance.value,
            "initial_magnetic_moments": initial_magnetic_moments,
            "pw": {
                "parameters": {
                    "SYSTEM": {
                        "degauss": self.degauss.value,
                        "smearing": self.smearing.value,
                        "tot_charge": self.tot_charge.value,
                    }
                },
            },
        }
        return parameters

    def set_panel_value(self, parameters):
        """Set a dictionary to set the value of the widgets in the panel.

        :param parameters: a dictionary of the values of all the widgets in the panel.
        """
        # TODO we shoudl first check the overide box, and then set the value
        # but discuss with others if we really need the orveride box
        pseudo_family = parameters.get("pseudo_family")
        self.pseudo_family_selector.value = pseudo_family
        # TODO support pseudodojo
        family, _version, functional, protocol = pseudo_family.split("/")
        self.pseudo_family_selector.dft_functional.value = functional
        self.pseudo_family_selector.protocol_selection.value = f"{family} {protocol}"
        self.kpoints_distance.value = parameters.get("kpoints_distance", 0.15)
        if parameters.get("pw") is not None:
            self.degauss.value = parameters["pw"]["parameters"]["SYSTEM"]["degauss"]
            self.smearing.value = parameters["pw"]["parameters"]["SYSTEM"]["smearing"]
            self.tot_charge.value = parameters["pw"]["parameters"]["SYSTEM"].get(
                "tot_charge", 0
            )
        self.magnetization._set_magnetization_values(
            parameters.get("initial_magnetic_moments", 0.0)
        )

    def reset(self):
        self.set_panel_value(DEFAULT_PARAMETERS)

    def _update_state(self):
        """Update the state of the panel."""
        if self.parent is not None:
            self.magnetization._update_widget(
                {"new": self.parent.parent.structure_step.confirmed_structure}
            )


class MagnetizationSettings(ipw.VBox):
    """Widget to set the initial magnetic moments for each kind names defined in the StructureData (StructureDtaa.get_kind_names())
    Usually these are the names of the elements in the StructureData
    (For example 'C' , 'N' , 'Fe' . However the StructureData can have defined kinds like 'Fe1' and 'Fe2')

    The widget generate a dictionary that can be used to set initial_magnetic_moments in the builder of PwBaseWorkChain

    Attributes:
        input_structure(StructureData): trait that containes the input_strucgure (confirmed structure from previous step)
    """

    def __init__(self, **kwargs):
        self.input_structure = StructureData()
        self.input_structure_labels = []
        self.description = ipw.HTML(
            "<b>Magnetization</b>: Input structure not confirmed"
        )
        self.kinds = self.create_kinds_widget()
        self.kinds_widget_out = ipw.Output()
        self.override = ipw.Checkbox(
            description="Override",
            indent=False,
            value=False,
        )
        super().__init__(
            children=[
                ipw.HBox(
                    [
                        # self.override,
                        self.description,
                        self.kinds_widget_out,
                    ],
                ),
            ],
            layout=ipw.Layout(justify_content="space-between"),
            **kwargs,
        )
        self.display_kinds()
        self.override.observe(self._disable_kinds_widgets, "value")

    def _disable_kinds_widgets(self, _=None):
        for i in range(len(self.kinds.children)):
            self.kinds.children[i].disabled = not self.override.value

    def reset(self):
        self.override.value = False
        if hasattr(self.kinds, "children") and self.kinds.children:
            for i in range(len(self.kinds.children)):
                self.kinds.children[i].value = 0.0

    def create_kinds_widget(self):
        if self.input_structure_labels:
            widgets_list = []
            for kind_label in self.input_structure_labels:
                kind_widget = ipw.BoundedFloatText(
                    description=kind_label,
                    min=-4,
                    max=4,
                    step=0.1,
                    value=0.0,
                    disabled=True,
                )
                widgets_list.append(kind_widget)
            kinds_widget = ipw.VBox(widgets_list)
        else:
            kinds_widget = ipw.VBox([])

        return kinds_widget

    def update_kinds_widget(self):
        self.input_structure_labels = self.input_structure.get_kind_names()
        self.kinds = self.create_kinds_widget()
        self.description.value = "<b> Magnetization</b>: "
        self.display_kinds()

    def display_kinds(self):
        import os

        if "PYTEST_CURRENT_TEST" not in os.environ and self.kinds:
            with self.kinds_widget_out:
                clear_output()
                display(self.kinds)

    def _update_widget(self, change):
        self.input_structure = change["new"]
        self.update_kinds_widget()

    def get_magnetization(self):
        """Method to generate the dictionary with the initial magnetic moments"""
        magnetization = {}
        for i in range(len(self.kinds.children)):
            magnetization[self.input_structure_labels[i]] = self.kinds.children[i].value
        return magnetization

    def _set_magnetization_values(self, initial_magnetic_moments):
        """Update used for conftest setting all magnetization to a value"""
        with self.hold_trait_notifications():
            if isinstance(initial_magnetic_moments, dict):
                for i in range(len(self.kinds.children)):
                    self.kinds.children[i].value = initial_magnetic_moments.get(
                        self.input_structure_labels[i], 0.0
                    )
            elif isinstance(initial_magnetic_moments, (int, float)):
                for i in range(len(self.kinds.children)):
                    self.kinds.children[i].value = initial_magnetic_moments
