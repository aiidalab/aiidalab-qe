"""Panel for Pdos plugin."""

import ipywidgets as ipw

from aiidalab_qe.common.panel import SettingsPanel

from .model import PdosModel

RYDBERG_TO_EV = 13.605703976


class Setting(SettingsPanel[PdosModel]):
    title = "PDOS"
    identifier = "pdos"

    def __init__(self, model: PdosModel, **kwargs):
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

        self.children = [
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 0px">
                    <h4>Settings</h4>
                </div>
            """),
            ipw.HTML("""
                <div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
                    By default, the tetrahedron method is used for PDOS calculation.
                    If required you can apply Gaussian broadening with a custom degauss
                    value.
                    <br>
                    For molecules and systems with localized orbitals, it is
                    recommended to use a custom degauss value.
                </div>
            """),
            ipw.HBox(
                children=[
                    self.nscf_kpoints_distance,
                    self.mesh_grid,
                ]
            ),
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
