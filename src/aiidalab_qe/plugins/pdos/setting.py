"""Panel for Pdos plugin."""

import ipywidgets as ipw

from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.panel import ConfigurationSettingsPanel

from .model import PdosConfigurationSettingsModel

RYDBERG_TO_EV = 13.605703976


class PdosConfigurationSettingPanel(
    ConfigurationSettingsPanel[PdosConfigurationSettingsModel],
):
    def __init__(self, model: PdosConfigurationSettingsModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_protocol_change,
            "protocol",
        )
        self._model.observe(
            self._on_kpoints_distance_change,
            "kpoints_distance",
        )

    def render(self):
        if self.rendered:
            return

        self.nscf_kpoints_distance = ipw.BoundedFloatText(
            min=0.001,
            step=0.01,
            description="NSCF K-points distance (1/Å):",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "kpoints_distance"),
            (self.nscf_kpoints_distance, "value"),
        )
        ipw.dlink(
            (self._model, "input_structure"),
            (self.nscf_kpoints_distance, "disabled"),
            lambda _: not self._model.has_pbc,
        )

        self.mesh_grid = ipw.HTML()
        ipw.dlink(
            (self._model, "mesh_grid"),
            (self.mesh_grid, "value"),
        )

        self.use_pdos_degauss = ipw.Checkbox(
            indent=False,
            description="Use custom PDOS degauss",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "use_pdos_degauss"),
            (self.use_pdos_degauss, "value"),
        )
        ipw.dlink(
            (self._model, "input_structure"),
            (self.use_pdos_degauss, "disabled"),
            lambda _: not self._model.has_pbc,
        )

        self.pdos_degauss = ipw.BoundedFloatText(
            min=0.00001,
            step=0.001,
            description="PDOS degauss (Ry):",
            disabled=True,
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "pdos_degauss"),
            (self.pdos_degauss, "value"),
        )
        ipw.dlink(
            (self._model, "use_pdos_degauss"),
            (self.pdos_degauss, "disabled"),
            lambda use_pdos_degauss: not use_pdos_degauss,
        )

        self.pdos_degauss_eV = ipw.HTML()
        ipw.dlink(
            (self._model, "pdos_degauss"),
            (self.pdos_degauss_eV, "value"),
            lambda degauss: f"({degauss * RYDBERG_TO_EV:.4f} eV)",
        )

        self.energy_grid_step = ipw.BoundedFloatText(
            min=0.001,
            step=0.001,
            description="Energy grid step (eV):",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "energy_grid_step"),
            (self.energy_grid_step, "value"),
        )

        self.children = [
            InAppGuide(identifier="pdos-settings"),
            ipw.HTML("""
                <div style="line-height: 140%; margin-bottom: 10px;">
                    <p style="margin-bottom: 10px;">
                        By default, the <b>tetrahedron method</b> is used for partial
                        density of states (PDOS) calculation. However, if you need more
                        control over the broadening, you can apply <b>Gaussian broadening</b>
                        by specifying a custom <b>degauss</b> value.
                    </p>
                    <p>
                        For systems involving <b>molecules</b> or <b>localized orbitals</b>,
                        it is recommended to use a <b>custom degauss value</b>. This
                        will provide a more accurate representation of the PDOS,
                        especially when the electronic states are localized.
                    </p>
                </div>
            """),
            ipw.HBox(
                children=[
                    self.nscf_kpoints_distance,
                    self.mesh_grid,
                ]
            ),
            self.energy_grid_step,
            self.use_pdos_degauss,
            ipw.HBox(
                children=[
                    self.pdos_degauss,
                    self.pdos_degauss_eV,
                ]
            ),
        ]

        self.rendered = True

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_protocol_change(self, _):
        self.refresh(specific="protocol")

    def _on_kpoints_distance_change(self, _):
        self.refresh(specific="mesh")
