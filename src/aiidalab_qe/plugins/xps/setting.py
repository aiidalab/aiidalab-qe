"""Panel for XPS plugin."""

import ipywidgets as ipw

from aiidalab_qe.app.configuration.model import ConfigurationModel
from aiidalab_qe.common.panel import SettingsPanel

from .model import BASE_URL


class Setting(SettingsPanel):
    title = "XPS Settings"
    identifier = "xps"

    def __init__(self, config_model: ConfigurationModel, **kwargs):
        super().__init__(config_model, **kwargs)
        ipw.dlink(
            (self._config_model, "input_structure"),
            (self._model, "input_structure"),
        )
        self._config_model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_pseudo_group_change,
            "pseudo_group",
        )

    def render(self):
        if self.rendered:
            return

        self.core_hole_treatment = ipw.ToggleButtons(
            options=[
                ("XCH(smear)", "xch_smear"),
                ("XCH(fixed)", "xch_fixed"),
                ("Full", "full"),
            ],
        )
        ipw.link(
            (self._model, "core_hole_treatment"),
            (self.core_hole_treatment, "value"),
        )

        self.pseudo_group = ipw.Dropdown(
            options=["pseudo_demo_pbe", "pseudo_demo_pbesol"],
            description="Group:",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "pseudo_group"),
            (self.pseudo_group, "value"),
        )

        self.core_levels_widget = ipw.VBox()

        self.structure_type = ipw.ToggleButtons(
            options=[
                ("Molecule", "molecule"),
                ("Crystal", "crystal"),
            ],
        )
        ipw.link(
            (self._model, "structure_type"),
            (self.structure_type, "value"),
        )

        self.supercell_min_parameter = ipw.FloatText(
            description="The minimum cell length (Å):",
            disabled=False,
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "supercell_min_parameter"),
            (self.supercell_min_parameter, "value"),
        )

        self.calc_binding_energy = ipw.Checkbox(
            description="Calculate binding energy: ",
            indent=False,
        )
        ipw.link(
            (self._model, "calc_binding_energy"),
            (self.calc_binding_energy, "value"),
        )

        self.children = [
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Structure</h4>
                </div>
            """),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px">
                    Below you can indicate if the material should be treated as a
                    molecule or a crystal.
                </div>"""),
            ipw.HBox([self.structure_type]),
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Core-Hole pseudopotential group</h4>
                </div>"""),
            ipw.HTML(f"""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px">
                    Please select a pseudopotential group, which provide the
                    ground-state and excited-state pseudopotentials for the element.
                    The pseudopotentials are downloaded from this <a href="{BASE_URL}">
                    repository</a>.
                </div>"""),
            self.pseudo_group,
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Select core-level</h4>
                </div>"""),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 6px; padding-bottom: 0px">
                    The list of core-levels to be considered for analysis.
                </div>"""),
            ipw.HBox([self.core_levels_widget]),
        ]

        self.rendered = True

    def update(self):
        if not self.updated:
            self._update()
            self.updated = True

    def reset(self):
        if not self._config_model.input_structure:
            self._unsubscribe()
        self._model.reset()
        self.updated = False

    def _on_input_structure_change(self, _=None):
        self._update(rebuild=True)

    def _on_pseudo_group_change(self, _=None):
        self._update(rebuild=True)

    def _update(self, rebuild=False):
        self._unsubscribe()
        self._show_loading()
        self._model.update()
        self._build_core_levels_widget(rebuild=rebuild)

    def _show_loading(self):
        if self.rendered:
            self.core_levels_widget.children = [self.loading_message]

    def _build_core_levels_widget(self, rebuild=False):
        if (
            not self.rendered
            or len(self.core_levels_widget.children) > 1
            and not rebuild
        ):
            return

        if not self._model.include or self._model.input_structure is None:
            return

        children = []

        kind_list = [Kind.symbol for Kind in self._model.input_structure.kinds]

        supported_core_levels = self._model.get_supported_core_levels()

        for element in kind_list:
            if element in supported_core_levels:
                for orbital in supported_core_levels[element]:
                    checkbox = ipw.Checkbox(
                        description=orbital,
                        indent=False,
                        layout=ipw.Layout(max_width="100%"),
                    )
                    link = ipw.link(
                        (self._model, "core_levels"),
                        (checkbox, "value"),
                        [
                            lambda cl, orbital=orbital: cl.get(orbital, False),
                            lambda v, orbital=orbital: {
                                **self._model.core_levels,
                                orbital: v,
                            },
                        ],
                    )
                    self.links.append(link)
                    children.append(checkbox)
            else:
                checkbox = ipw.Checkbox(
                    description=f"{element}, not supported by the selected pseudo group",
                    indent=False,
                    disabled=True,
                    value=False,
                    style={"description_width": "initial"},
                    layout=ipw.Layout(max_width="100%"),
                )
                children.append(checkbox)

        self.core_levels_widget.children = children
