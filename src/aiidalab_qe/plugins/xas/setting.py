"""Panel for XAS plugin."""

import ipywidgets as ipw

from aiidalab_qe.common.panel import SettingsPanel

from .model import XasModel


class XasSettings(SettingsPanel[XasModel]):
    title = "XAS"
    identifier = "xas"

    # TODO: The element selection should lock the "Confirm" button if no elements have been selected for XAS calculation.

    def __init__(self, model: XasModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )

    def render(self):
        if self.rendered:
            return

        self.core_hole_treatments_widget = ipw.VBox(layout=ipw.Layout(width="100%"))

        # self.structure_type = ipw.ToggleButtons()
        # ipw.dlink(
        #     (self._model, "structure_type_options"),
        #     (self.structure_type, "options"),
        # )
        # ipw.link(
        #     (self._model, "structure_type"),
        #     (self.structure_type, "value"),
        # )

        self.supercell_min_parameter = ipw.FloatText(
            description="The minimum cell length (Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "supercell_min_parameter"),
            (self.supercell_min_parameter, "value"),
        )

        self.children = [
            # I will leave these objects here for now (15/11/23), but since the
            # calculation of molecular systems is not really supported (neither in
            # terms of XAS nor the main App itself) we should not present this option
            # that essentially does nothing.
            # ipw.HTML("""
            #     <div style="padding-top: 0px; padding-bottom: 0px">
            #         <h4>Structure</h4>
            #     </div>
            # """),
            # ipw.HTML("""
            #     <div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
            #         Below you can indicate if the material should be treated as a
            #         molecule or a crystal.
            #     </div>
            # """),
            # ipw.HBox(
            #     children=[self.structure_type],
            # ),
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Element and Core-Hole Treatment Setting.</h4>
                </div>
            """),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
                    To select elements for calculation of K-edge spectra:
                    <br>
                    (1) Tick the checkbox for each element symbol to select the element
                    for calculation.<br>
                    (2) Select the core-hole treatment scheme from the dropdown box.
                    <br>
                    <br>
                    There are three supported options for core-hole treatment:
                    <br>
                    - FCH: Remove one electron from the system (any occupations scheme).
                    <br>
                    - XCH (Smearing): places the excited electron into the conduction
                    band (smeared occupations).
                    <br>
                    - XCH (Fixed): places the excited electron into the conduction band
                    (fixed occupations).
                    <br>
                    <br>
                    For XAS calculations of most elements, the FCH treatment is
                    recommended, however in some cases the XCH treatment should be used instead.
                    <br>
                    The recommended setting will be shown for each available element.
                    Note that only elements for which core-hole pseudopotential sets
                    are available will be shown.
                    <br>
                </div>
            """),
            ipw.HBox(
                children=[self.core_hole_treatments_widget],
                layout=ipw.Layout(width="95%"),
            ),
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Cell size</h4>
                </div>
            """),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
                    Define the minimum cell length in angstrom for the resulting
                    supercell, and thus all output structures. The default value of 8.0
                    angstrom will be used if no input is given. Setting this value to
                    0.0 will instruct the CF to not scale up the input structure.
                </div>
            """),
            ipw.HBox(
                children=[self.supercell_min_parameter],
            ),
        ]

        self.rendered = True

        self.refresh(specific="widgets")

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def update(self, specific=""):
        if self.updated:
            return
        self._show_loading()
        if not self._model.loaded_from_process or specific and specific != "widgets":
            self._model.update(specific)
        self._build_core_hole_treatments_widget()
        self.updated = True

    def _show_loading(self):
        if self.rendered:
            self.core_hole_treatments_widget.children = [self.loading_message]

    def _build_core_hole_treatments_widget(self):
        if not self.rendered:
            return

        children = []

        info = "Recommended: <b>{recommended}</b> (PBE Core-Hole Pseudopotential)"

        for kind_name in self._model.get_kind_names():
            checkbox = ipw.Checkbox(
                description=kind_name,
                disabled=False,
                style={"description_width": "initial"},
                layout=ipw.Layout(width="7%"),
            )
            link = ipw.link(
                (self._model, "kind_names"),
                (checkbox, "value"),
                [
                    lambda kinds, kind_name=kind_name: kinds.get(kind_name, False),
                    lambda value, kind_name=kind_name: {
                        **self._model.kind_names,
                        kind_name: value,
                    },
                ],
            )
            self.links.append(link)

            treatment_selector = ipw.Dropdown(layout=ipw.Layout(width="15%"))
            options_link = ipw.dlink(
                (self._model, "core_hole_treatments_options"),
                (treatment_selector, "options"),
            )
            value_link = ipw.link(
                (self._model, "core_hole_treatments"),
                (treatment_selector, "value"),
                [
                    lambda treatments, kind_name=kind_name: treatments.get(
                        kind_name, "full"
                    ),
                    lambda value, kind_name=kind_name: {
                        **self._model.core_hole_treatments,
                        kind_name: value,
                    },
                ],
            )
            self.links.extend([options_link, value_link])

            widget = ipw.HBox(
                children=[
                    checkbox,
                    treatment_selector,
                    ipw.HTML(
                        value=info.format(
                            recommended=self._model.get_recommendation(kind_name)
                        ),
                        layout=ipw.Layout(width="78%"),
                    ),
                ],
                layout=ipw.Layout(width="100%"),
            )

            children.append(widget)

        if not children:
            children = [ipw.HTML("<b>No elements available for XAS calculations</b>")]

        self.core_hole_treatments_widget.children = children

        # For reference:
        # This is the whole widget:
        # print(f"{self.core_hole_treatments_widget}\n")

        # This is the tuple of selected element and core-hole treatment:
        # print(f"{self.core_hole_treatments_widget.children[0]}\n")

        # This is the checkbox for the element, giving element name and whether to add it to the elements list
        # print(f"{self.core_hole_treatments_widget.children[0].children[0]}\n")
        # print(f"{self.core_hole_treatments_widget.children[0].children[0].value}\n")
        # print(f"{self.core_hole_treatments_widget.children[0].children[0].description}\n")

        # This is the dropdown for the core-hole treatment option:
        # print(f"{self.core_hole_treatments_widget.children[0].children[1]}\n")
        # print(f"{self.core_hole_treatments_widget.children[0].children[1].value}\n")
