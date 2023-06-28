import ipywidgets as ipw
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain

from aiidalab_qe.app.panel import Panel
from aiidalab_qe.app.pseudos import PseudoFamilySelector


class AdvancedSettings(Panel):
    title = "Advanced Settings"

    subtitle = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Advanced Settings</h4></div>"""
    )
    subdescription = ipw.HTML(
        """Select the advanced settings for the <b>pw.x</b> code."""
    )

    smearing_description = ipw.HTML(
        """<p>
        The smearing type and width is set by the chosen <b>protocol</b>.
        Tick the box to override the default, not advised unless you've mastered <b>smearing effects</b> (click <a href="http://theossrv1.epfl.ch/Main/ElectronicTemperature"
        target="_blank">here</a> for a discussion).
    </p>"""
    )
    kpoints_distance_description = ipw.HTML(
        """<div>
        The k-points mesh density of the SCF calculation is set by the <b>protocol</b>.
        The value below represents the maximum distance between the k-points in each direction of reciprocal space.
        Tick the box to override the default, smaller is more accurate and costly. </div>"""
    )
    tot_charge_description = ipw.HTML("""<p></p>""")

    def __init__(self, **kwargs):
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
            description="Smearing type:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.degauss = ipw.FloatText(
            value=0.01,
            step=0.005,
            description="Smearing width (Ry):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.kpoints_distance = ipw.FloatText(
            value=0.15,
            step=0.05,
            description="K-points distance (1/Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.tot_charge = ipw.BoundedFloatText(
            value=0,
            min=-3,
            max=3,
            step=0.01,
            disabled=False,
            description="Total charge:",
            style={"description_width": "initial"},
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
            self.subdescription,
            # self.tot_charge_description,
            self.tot_charge,
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

        parameters = {
            "pseudo_family": self.pseudo_family_selector.value,
            "kpoints_distance": self.kpoints_distance.value,
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
        self.pseudo_family_selector.override_protocol_pseudo_family.value = True
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
