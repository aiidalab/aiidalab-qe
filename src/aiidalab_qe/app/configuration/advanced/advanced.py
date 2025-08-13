"""Widgets for the submission of bands work chains.

Authors: AiiDAlab team
"""

import ipywidgets as ipw

from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ConfigurationSettingsPanel

from .convergence import (
    ConvergenceConfigurationSettingsModel,
    ConvergenceConfigurationSettingsPanel,
)
from .hubbard import (
    HubbardConfigurationSettingsModel,
    HubbardConfigurationSettingsPanel,
)
from .magnetization import (
    MagnetizationConfigurationSettingsModel,
    MagnetizationConfigurationSettingsPanel,
)
from .model import AdvancedConfigurationSettingsModel
from .pseudos import (
    PseudosConfigurationSettingsModel,
    PseudosConfigurationSettingsPanel,
)
from .smearing import (
    SmearingConfigurationSettingsModel,
    SmearingConfigurationSettingsPanel,
)


class AdvancedConfigurationSettingsPanel(
    ConfigurationSettingsPanel[AdvancedConfigurationSettingsModel],
):
    def __init__(self, model: AdvancedConfigurationSettingsModel, **kwargs):
        super().__init__(
            model=model,
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            **kwargs,
        )

        self._model.observe(
            self._on_spin_type_change,
            "spin_type",
        )

        convergence_model = ConvergenceConfigurationSettingsModel()
        self.convergence = ConvergenceConfigurationSettingsPanel(
            model=convergence_model
        )
        model.add_model("convergence", convergence_model)

        smearing_model = SmearingConfigurationSettingsModel()
        self.smearing = SmearingConfigurationSettingsPanel(model=smearing_model)
        model.add_model("smearing", smearing_model)

        pseudos_model = PseudosConfigurationSettingsModel()
        self.pseudos = PseudosConfigurationSettingsPanel(model=pseudos_model)
        model.add_model("pseudos", pseudos_model)

        magnetization_model = MagnetizationConfigurationSettingsModel()
        self.magnetization = MagnetizationConfigurationSettingsPanel(
            model=magnetization_model,
        )
        model.add_model("magnetization", magnetization_model)

        hubbard_model = HubbardConfigurationSettingsModel()
        self.hubbard = HubbardConfigurationSettingsPanel(model=hubbard_model)
        model.add_model("hubbard", hubbard_model)

        self.sub_settings = {
            "convergence": self.convergence,
            "smearing": self.smearing,
            "magnetization": self.magnetization,
            "hubbard": self.hubbard,
            "pseudos": self.pseudos,
        }

    def render(self):
        if self.rendered:
            return

        self.reset_to_defaults_button = ipw.Button(
            description="Reset to defaults",
            button_style="primary",
            icon="undo",
            layout=ipw.Layout(width="fit-content"),
        )
        self.reset_to_defaults_button.on_click(self._on_reset_to_defaults_button_click)

        self.clean_workdir = ipw.Checkbox(
            description="Delete the work directory after the calculation",
            indent=False,
            layout=ipw.Layout(width="fit-content", margin="5px 2px"),
        )
        ipw.link(
            (self._model, "clean_workdir"),
            (self.clean_workdir, "value"),
        )

        # Total change setting
        self.total_charge = ipw.BoundedFloatText(
            min=-3,
            max=3,
            step=0.01,
            description="Total charge:",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "total_charge"),
            (self.total_charge, "value"),
        )

        # Van der Waals setting
        self.van_der_waals = ipw.Dropdown(
            description="Van der Waals correction:",
            style={"description_width": "150px"},
        )
        ipw.dlink(
            (self._model, "van_der_waals_options"),
            (self.van_der_waals, "options"),
        )
        ipw.link(
            (self._model, "van_der_waals"),
            (self.van_der_waals, "value"),
        )

        self.advanced_tabs = ipw.Tab()
        self.advanced_tabs.observe(
            self._on_advanced_tab_change,
            "selected_index",
        )

        self.children = [
            InAppGuide(identifier="advanced-settings"),
            self.reset_to_defaults_button,
            self.clean_workdir,
            self.total_charge,
            self.van_der_waals,
            ipw.HTML("<hr>"),
            self.advanced_tabs,
        ]

        self.rendered = True

        self.refresh()

        self._update_tabs()

    def _on_spin_type_change(self, _):
        self._update_tabs()

    def _update_tabs(self):
        if not self.rendered:
            return
        children = []
        titles = []
        for identifier, model in self._model.get_models():
            subsetting = self.sub_settings[identifier]
            if identifier == "convergence":
                subsetting.render()
            if identifier == "magnetization":
                if self._model.spin_type != "collinear":
                    continue
                subsetting.render()
            titles.append(model.title)
            children.append(subsetting)
        self.advanced_tabs.children = children
        for i, title in enumerate(titles):
            self.advanced_tabs.set_title(i, title)

    def _on_advanced_tab_change(self, change):
        self.advanced_tabs.children[change["new"]].render()

    def _on_reset_to_defaults_button_click(self, _):
        self._reset()

    def _reset(self):
        self._model.reset()
        for _, model in self._model.get_models():
            model.reset()
