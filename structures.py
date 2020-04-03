"""Widgets for the upload and selection of structure data.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
import os
import sys
import tempfile

import aiida
import ase
import nglview
import traitlets
import ipywidgets as ipw

from wizard import WizardApp, WizardAppStep


class StructureSelectionStep(ipw.VBox, WizardAppStep):
    """Integrated widget for the selection of structures from different sources."""

    structure = traitlets.Instance(aiida.orm.StructureData, allow_none=True)
    confirmed_structure = traitlets.Instance(aiida.orm.StructureData, allow_none=True)

    def __init__(self, manager, description=None, **kwargs):
        self.manager = manager

        if description is None:
            description = ipw.Label(
                'Select a structure from one of the following sources and then '
                'click "Confirm" to go to the next step.')
        self.description = description

        self.structure_name_text = ipw.Text(
            placeholder='[No structure selected]',
            description='Selected:',
            disabled=True,
            layout=ipw.Layout(width='auto', flex="1 1 auto"),
        )

        self.confirm_button = ipw.Button(
            description='Confirm',
            tooltip="Confirm the currently selected structure and go to the next step.",
            button_style='success',
            icon='check-circle',
            disabled=True,
            layout=ipw.Layout(width='auto'),
        )
        self.confirm_button.on_click(self.confirm)

        # Create directional link from the (read-only) 'structure_node' traitlet of the
        # structure manager to our 'structure' traitlet:
        ipw.dlink((manager, 'structure_node'), (self, 'structure'))

        super().__init__(
            children=[self.description, self.manager, self.structure_name_text, self.confirm_button], **kwargs)

    @traitlets.default('state')
    def _default_state(self):
        return WizardApp.State.READY

    def _update_state(self):
        if self.structure is None:
            if self.confirmed_structure is None:
                self.state = WizardApp.State.READY
            else:
                self.state = WizardApp.State.FAIL
        else:
            if self.confirmed_structure is None:
                self.state = WizardApp.State.CONFIGURED
            else:
                self.state = WizardApp.State.SUCCESS

    @traitlets.observe('structure')
    def _observe_structure(self, change):
        structure = change['new']
        with self.hold_trait_notifications():
            if structure is None:
                self.structure_name_text.value = ""
            else:
                self.structure_name_text.value = str(self.structure.get_formula())
            self._update_state()

    @traitlets.observe('confirmed_structure')
    def _observe_confirmed_structure(self, _):
        with self.hold_trait_notifications():
            self._update_state()

    @traitlets.observe('state')
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            state = change['new']
            self.confirm_button.disabled = state != WizardApp.State.CONFIGURED
            self.manager.disabled = state is WizardApp.State.SUCCESS

    def confirm(self, _=None):
        self.confirmed_structure = self.structure

    def reset(self):  # unconfirm
        self.confirmed_structure = None
