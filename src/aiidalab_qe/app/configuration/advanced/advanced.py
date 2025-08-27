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
from .general import (
    GeneralConfigurationSettingsModel,
    GeneralConfigurationSettingsPanel,
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
            layout={
                "justify_content": "space-between",
                "grid_gap": "4px",
                **kwargs.get("layout", {}),
            },
            **kwargs,
        )

        general_model = GeneralConfigurationSettingsModel()
        self.general = GeneralConfigurationSettingsPanel(model=general_model)
        model.add_model("general", general_model)

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

        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_spin_type_change,
            "spin_type",
        )

        self.sub_settings: dict[str, ConfigurationSettingsPanel] = {
            "general": self.general,
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

        self.advanced_tabs = ipw.Tab()
        self.advanced_tabs.observe(
            self._on_advanced_tab_change,
            "selected_index",
        )

        self.children = [
            InAppGuide(identifier="advanced-settings"),
            self.reset_to_defaults_button,
            self.advanced_tabs,
        ]

        self.rendered = True

        self._update_tabs()

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_spin_type_change(self, _):
        self._update_tabs()

    def _update_tabs(self):
        if not self.rendered:
            return
        children = []
        titles = []
        for identifier, model in self._model.get_models():
            subsetting = self.sub_settings[identifier]
            if identifier == "general":
                subsetting.render()
            if not model.include:
                continue
            titles.append(model.title)
            children.append(subsetting)
        self.advanced_tabs.children = children
        for i, title in enumerate(titles):
            self.advanced_tabs.set_title(i, title)
        self.advanced_tabs.selected_index = 0

    def _on_advanced_tab_change(self, change):
        tab: ConfigurationSettingsPanel = self.advanced_tabs.children[change["new"]]
        tab.render()

    def _on_reset_to_defaults_button_click(self, _):
        self._reset()

    def _reset(self):
        self._model.reset()
        for _, model in self._model.get_models():
            model.reset()
