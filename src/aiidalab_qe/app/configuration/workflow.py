"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import typing as t

import ipywidgets as ipw

from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.panel import SettingsModel, SettingsPanel


class WorkChainSettings(SettingsPanel):
    title = "Basic Settings"
    identifier = "workchain"

    def fetch_setting_entries(
        self,
        register_setting_callback: t.Callable[[str, SettingsPanel], None],
        update_tabs_callback: t.Callable[[], None],
    ):
        self.property_children = [ipw.HTML("Select which properties to calculate:")]

        outlines = get_entry_items("aiidalab_qe.properties", "outline")
        models = get_entry_items("aiidalab_qe.properties", "model")
        settings = get_entry_items("aiidalab_qe.properties", "setting")
        for identifier in settings:
            model: SettingsModel = models[identifier]()
            self._config_model.add_model(identifier, model)

            outline = outlines[identifier]()
            info = ipw.HTML()
            ipw.link(
                (model, "include"),
                (outline.include, "value"),
            )

            if identifier == "bands":
                ipw.dlink(
                    (self._config_model, "has_pbc"),
                    (outline.include, "disabled"),
                    lambda periodic: not periodic,
                )

            def toggle_plugin(change, identifier=identifier, info=info):
                if change["new"]:
                    info.value = (
                        f"Customize {identifier} settings in the panel above if needed"
                    )
                else:
                    info.value = ""
                update_tabs_callback()

            model.observe(
                toggle_plugin,
                "include",
            )

            self.property_children.append(
                ipw.HBox(
                    children=[
                        outline,
                        info,
                    ]
                )
            )

            register_setting_callback(identifier, settings[identifier])

    def update(self):
        if self.updated:
            return
        self._model.update()
        self.updated = True

    def render(self):
        if self.rendered:
            return

        # RelaxType: degrees of freedom in geometry optimization
        self.relax_type_help = ipw.HTML()
        ipw.dlink(
            (self._model, "relax_type_help"),
            (self.relax_type_help, "value"),
        )
        self.relax_type = ipw.ToggleButtons(
            options=[
                ("Structure as is", "none"),
                ("Atomic positions", "positions"),
                ("Full geometry", "positions_cell"),
            ],
        )
        ipw.dlink(
            (self._model, "relax_type_options"),
            (self.relax_type, "options"),
        )
        ipw.link(
            (self._model, "relax_type"),
            (self.relax_type, "value"),
        )

        # SpinType: magnetic properties of material
        self.spin_type = ipw.ToggleButtons(
            options=[("Off", "none"), ("On", "collinear")],
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "spin_type"),
            (self.spin_type, "value"),
        )

        # ElectronicType: electronic properties of material
        self.electronic_type = ipw.ToggleButtons(
            options=[("Metal", "metal"), ("Insulator", "insulator")],
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "electronic_type"),
            (self.electronic_type, "value"),
        )

        # Work chain protocol
        self.protocol = ipw.ToggleButtons(
            options=["fast", "moderate", "precise"],
        )
        ipw.link(
            (self._model, "protocol"),
            (self.protocol, "value"),
        )

        self.children = [
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Structure</h4>
                </div>
            """),
            self.relax_type_help,
            self.relax_type,
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
                    Below you can indicate both if the material should be treated as an
                    insulator or a metal (if in doubt, choose "Metal"), and if it
                    should be studied with magnetization/spin polarization, switch
                    magnetism On or Off (On is at least twice more costly).
                </div>
            """),
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Electronic Type:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.electronic_type,
                ]
            ),
            ipw.HBox(
                children=[
                    ipw.Label(
                        "Magnetism:",
                        layout=ipw.Layout(justify_content="flex-start", width="120px"),
                    ),
                    self.spin_type,
                ]
            ),
            *self.property_children,
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Protocol</h4>
                </div>
            """),
            ipw.HTML("Select the protocol:", layout=ipw.Layout(flex="1 1 auto")),
            self.protocol,
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 6px; padding-bottom: 0px">
                    The "moderate" protocol represents a trade-off between accuracy and
                    speed. Choose the "fast" protocol for a faster calculation with
                    less precision and the "precise" protocol to aim at best accuracy
                    (at the price of longer/costlier calculations).
                </div>
            """),
        ]

        self.rendered = True

    def reset(self):
        self._model.reset()
        self.updated = False
