"""Widgets for the upload and selection of structure data.

Authors: AiiDAlab team
"""

import pathlib
import aiida
import ipywidgets as ipw
import traitlets
from aiidalab_widgets_base import WizardAppWidgetStep

# The Examples list of (name, file) tuple curretly passed to
# StructureExamplesWidget.
file_path = pathlib.Path(__file__).parent
Examples = [
    ("Silicon (diamond)", file_path / "examples" / "Si.xyz"),
    ("Silicon oxide", file_path / "examples" / "SiO2.xyz"),
    ("Diamond", file_path / "examples" / "diamond.cif"),
    ("Gallium arsenide", file_path / "examples" / "GaAs.xyz"),
    ("Gold (fcc)", file_path / "examples" / "Au.cif"),
    ("Cobalt (hcp)", file_path / "examples" / "Co.cif"),
]


class StructureSelectionStep(ipw.VBox, WizardAppWidgetStep):
    """Integrated widget for the selection of structures from different sources."""

    structure = traitlets.Instance(aiida.orm.StructureData, allow_none=True)
    confirmed_structure = traitlets.Instance(aiida.orm.StructureData, allow_none=True)

    def __init__(self, manager, description=None, **kwargs):
        self.manager = manager

        if description is None:
            description = ipw.HTML(
                """
                <p>Select a structure from one of the following sources and then click
                "Confirm" to go to the next step. </p><i class="fa fa-exclamation-circle"
                aria-hidden="true"></i> Currently only three-dimensional structures are
                supported.
                """
            )
        self.description = description

        self.structure_name_text = ipw.Text(
            placeholder="[No structure selected]",
            description="Selected:",
            disabled=True,
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
        )

        self.confirm_button = ipw.Button(
            description="Confirm",
            tooltip="Confirm the currently selected structure and go to the next step.",
            button_style="success",
            icon="check-circle",
            disabled=True,
            layout=ipw.Layout(width="auto"),
        )
        self.confirm_button.on_click(self.confirm)
        self.message_area = ipw.HTML()

        # Create directional link from the (read-only) 'structure_node' traitlet of the
        # structure manager to our 'structure' traitlet:
        ipw.dlink((manager, "structure_node"), (self, "structure"))

        super().__init__(
            children=[
                self.description,
                self.manager,
                self.structure_name_text,
                self.message_area,
                self.confirm_button,
            ],
            **kwargs
        )

    @traitlets.default("state")
    def _default_state(self):
        return self.State.INIT

    def _update_state(self):
        if self.structure is None:
            if self.confirmed_structure is None:
                self.state = self.State.READY
            else:
                self.state = self.State.SUCCESS
        else:
            if self.confirmed_structure is None:
                self.state = self.State.CONFIGURED
            else:
                self.state = self.State.SUCCESS

    @traitlets.observe("structure")
    def _observe_structure(self, change):
        structure = change["new"]
        with self.hold_trait_notifications():
            if structure is None:
                self.structure_name_text.value = ""
                self.message_area.value = ""
            else:
                self.structure_name_text.value = str(self.structure.get_formula())
            self._update_state()

    @traitlets.observe("confirmed_structure")
    def _observe_confirmed_structure(self, _):
        with self.hold_trait_notifications():
            self._update_state()

    @traitlets.observe("state")
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            state = change["new"]
            self.confirm_button.disabled = state != self.State.CONFIGURED
            self.manager.disabled = state is self.State.SUCCESS

    def confirm(self, _=None):
        self.manager.store_structure()
        self.confirmed_structure = self.structure
        self.message_area.value = ""

    def can_reset(self):
        return self.confirmed_structure is not None

    def reset(self):  # unconfirm
        self.confirmed_structure = None
