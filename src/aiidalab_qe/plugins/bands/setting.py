"""Panel for Bands plugin."""

import ipywidgets as ipw

from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ConfigurationSettingsPanel
from aiidalab_qe.plugins.bands.model import BandsConfigurationSettingsModel


class BandsConfigurationSettingsPanel(
    ConfigurationSettingsPanel[BandsConfigurationSettingsModel],
):
    def render(self):
        if self.rendered:
            return

        self.projwfc_bands = ipw.Checkbox(
            indent=False,
            description="Fat bands calculation",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "projwfc_bands"),
            (self.projwfc_bands, "value"),
        )

        self.children = [
            InAppGuide(identifier="bands-settings"),
            ipw.HTML("""
                <div style="line-height: 140%;">
                    <p style="margin-bottom: 10px;">
                        The band structure workflow will automatically detect the
                        default path in reciprocal space using the
                        <a
                            href="https://www.materialscloud.org/work/tools/seekpath"
                            target="_blank"
                        >SeeK-path tool</a>.
                    </p>
                    <p>
                        Fat Bands is a band structure plot that includes the angular
                        momentum contributions from specific atoms or orbitals to each
                        energy band. The thickness of the bands represents the strength
                        of these contributions, providing insight into the electronic
                        structure.
                    </p>
                </div>
            """),
            self.projwfc_bands,
        ]

        self.rendered = True
