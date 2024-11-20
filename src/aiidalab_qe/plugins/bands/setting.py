"""Panel for Bands plugin."""

import ipywidgets as ipw

from aiidalab_qe.common.panel import SettingsPanel
from aiidalab_qe.plugins.bands.model import BandsModel


class BandsSettings(SettingsPanel[BandsModel]):
    title = "Bands Structure"
    identifier = "bands"

    def render(self):
        if self.rendered:
            return

        self.projwfc_bands = ipw.Checkbox(
            description="Fat bands calculation",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "projwfc_bands"),
            (self.projwfc_bands, "value"),
        )

        self.children = [
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Settings</h4>
                </div>
            """),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px">
                    The band structure workflow will automatically detect the default
                    path in reciprocal space using the
                    <a href="https://www.materialscloud.org/work/tools/seekpath" target="_blank">SeeK-path tool</a>.
                    <br><br>
                    Fat Bands is a band structure plot that includes the angular
                    momentum contributions from specific atoms or orbitals to each
                    energy band. The thickness of the bands represents the strength
                    of these contributions, providing insight into the electronic
                    structure.
                </div>
            """),
            self.projwfc_bands,
        ]

        self.rendered = True
