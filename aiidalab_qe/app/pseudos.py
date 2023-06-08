# -*- coding: utf-8 -*-
import re

import ipywidgets as ipw
import traitlets

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS


class PseudoFamilySelector(ipw.VBox):
    title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Accuracy and precision</h4></div>"""
    )

    description = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        The exchange-correlation functional and pseudopotential library is set by
        the <b>protocol</b> configured in the "Workflow" tab. Here you can
        override the defaults if desired.</div>""",
        layout=ipw.Layout(max_width="60%"),
    )

    pseudo_family_prompt = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; opacity:0.5;">
        <b><a href="https://www.materialscloud.org/discover/sssp/table/efficiency"
        target="_blank">Pseudopotential family</b></div>"""
    )
    pseudo_family_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px; opacity:0.5;">
        If you are unsure, select 'SSSP efficiency', which for
        most calculations will produce sufficiently accurate results at
        comparatively small computational costs. If your calculations require a
        higher accuracy, select 'SSSP accuracy' or 'PseudoDojo stringent', which will be computationally
        more expensive. SSSP is the standard solid-state pseudopotentials.
        The PseudoDojo used here has the SR relativistic type.</div>"""
    )

    dft_functional_prompt = ipw.HTML(
        """
        <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; opacity:0.5;"><b>
        Exchange-correlation  functional</b></div>"""
    )
    dft_functional_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px; opacity:0.5;">
        The exchange-correlation energy is calculated using this functional. We currently provide support for two
        well-established generalised gradient approximation (GGA) functionals:
        PBE and PBEsol.</div>"""
    )
    disabled = traitlets.Bool()

    value = traitlets.Unicode(
        default_value=DEFAULT_PARAMETERS["pseudo_family"],
    )

    def __init__(self, **kwargs):
        # Enable manual setting of the pseudopotential family
        self.set_pseudo_family_prompt = ipw.HTML("<b>&nbsp;&nbsp;Override&nbsp;</b>")
        self.override_protocol_pseudo_family = ipw.Checkbox(
            description="",
            indent=False,
            value=False,
            layout=ipw.Layout(max_width="10%"),
        )
        self.set_pseudo_family_box = ipw.HBox(
            [self.set_pseudo_family_prompt, self.override_protocol_pseudo_family],
            layout=ipw.Layout(max_width="20%"),
        )
        self.show_ui = ipw.Valid(value=True)
        self.override_protocol_pseudo_family.observe(self.set_show_ui, "value")
        self.override_protocol_pseudo_family.observe(self.set_text_color, "value")
        self.override_protocol_pseudo_family.observe(self.set_value_trait, "value")

        # Choose the DFT functional
        self.dft_functional = ipw.Dropdown(
            options=["PBE", "PBEsol"],
            value=DEFAULT_PARAMETERS["pseudo_family"].split("/")[2],
            style={"description_width": "initial"},
        )
        self.dft_functional.observe(self.set_value_trait, "value")
        #
        pseudo_family_type = DEFAULT_PARAMETERS["pseudo_family"].split("/")[0]
        if pseudo_family_type.upper() == "SSSP":
            pseudo_family_type += (
                " " + DEFAULT_PARAMETERS["pseudo_family"].split("/")[-1]
            )
        elif pseudo_family_type.upper() == "PSEUDODOJO":
            pseudo_family_type = "PseudoDojo " + pseudo_family_type.split("_")[-2]
        self.protocol_selection = ipw.ToggleButtons(
            options=[
                "SSSP efficiency",
                "SSSP precision",
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ],
            value=pseudo_family_type,
            layout=ipw.Layout(max_width="100%"),
        )
        self.protocol_selection.observe(self.set_value_trait, "value")

        self.dft_functional_box = ipw.VBox(
            children=[
                self.dft_functional_prompt,
                self.dft_functional,
                self.dft_functional_help,
            ],
            layout=ipw.Layout(max_width="40%"),
        )
        self.pseudo_protocol_box = ipw.VBox(
            children=[
                self.pseudo_family_prompt,
                self.protocol_selection,
                self.pseudo_family_help,
            ],
            layout=ipw.Layout(max_width="60%"),
        )
        ipw.dlink((self.show_ui, "value"), (self.protocol_selection, "disabled"))
        ipw.dlink((self.show_ui, "value"), (self.dft_functional, "disabled"))

        self.set_value_trait()

        super().__init__(
            children=[
                self.title,
                ipw.HBox(
                    [self.description, self.set_pseudo_family_box],
                    layout=ipw.Layout(height="50px", justify_content="space-between"),
                ),
                ipw.HBox([self.dft_functional_box, self.pseudo_protocol_box]),
            ]
        )

    def set_value_trait(self, _=None):
        if self.override_protocol_pseudo_family.value:
            if "PseudoDojo" in self.protocol_selection.value:
                pseudo_family = f"PseudoDojo/0.4/{self.dft_functional.value}/SR/{self.protocol_selection.value.split()[1]}/upf"
            elif "SSSP" in self.protocol_selection.value:
                pseudo_family = f"SSSP/1.2/{self.dft_functional.value}/{self.protocol_selection.value.split()[1]}"
            else:
                raise ValueError(
                    f"Unknown pseudo family {self.override_protocol_pseudo_family.value}"
                )
            self.value = pseudo_family
        else:
            self.value = DEFAULT_PARAMETERS["pseudo_family"]

    def set_show_ui(self, change):
        self.show_ui.value = not change.new

    def set_text_color(self, change):
        opacity = 1.0 if change.new else 0.5

        for html in (
            self.pseudo_family_prompt,
            self.pseudo_family_help,
            self.dft_functional_help,
            self.dft_functional_prompt,
        ):
            old_opacity = re.match(
                r"[\s\S]+opacity:([\S]+);[\S\s]+", html.value
            ).groups()[0]
            html.value = html.value.replace(
                f"opacity:{old_opacity};", f"opacity:{opacity};"
            )

    def reset(self):
        self.protocol_selection.value = "SSSP efficiency"
