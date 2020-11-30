import ipywidgets as ipw
import traitlets

from aiida.orm import StructureData, Float, Str, BandsData
from aiida.plugins import WorkflowFactory
from aiida.engine import submit

from wizard import WizardApp
from codes import CodeSubmitWidget
from util import load_default_parameters


class ComputeBandsSubmitWidget(CodeSubmitWidget):

    input_structure = traitlets.Instance(StructureData, allow_none=True)
    band_structure = traitlets.Instance(BandsData, allow_none=True)

    def __init__(self, **kwargs):
        description = ipw.Label(
            'Specify the parameters and options for the calculation and then click on "Submit".')
        super().__init__(description=description, **kwargs)

    def _update_state(self):
        if self.process is None:
            if self.input_structure is None:
                if self.band_structure is None:
                    self.state = WizardApp.State.INIT
                else:
                    self.state = WizardApp.State.FAIL
            else:
                if self.band_structure is None:
                    self.state = WizardApp.State.READY
                else:
                    self.state = WizardApp.State.SUCCESS
        else:
            super()._update_state()

    @traitlets.observe('input_structure')
    def _observe_input_structure(self, change):
        self._update_state()

    @traitlets.observe('band_structure')
    def _observe_band_structure(self, change):
        self._update_state()

    @traitlets.validate('process')
    def _validate_process(self, proposal):
        process_node = proposal['value']
        builder = process_node.get_builder_restart()
        try:
            # Check that parameters are consistent with what we would expect for the app.
            # This ensures that both processes we create with the app and those that are
            # loaded are consistent.
            assert process_node.process_label == 'PwBandsWorkChain'
            assert builder.scf.pseudo_family.value == builder.bands.pseudo_family.value
            assert dict(load_default_parameters()) \
                    == dict(builder.scf.pw['parameters']) \
                    == dict(builder.bands.pw['parameters'])
            assert builder.scf.kpoints_distance == Float(0.8)
        except AssertionError as error:
            raise traitlets.TraitError(f"Unable to set process, parameters are inconsistent: {error}")

        return proposal['value']

    @traitlets.observe('process')
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            super()._observe_process(change)
            process_node = change['new']
            if process_node is None:
                self.band_structure = None
            else:
                # Restore the code parameters
                builder = process_node.get_builder_restart()
                self.pseudo_family_selection.value = builder.scf.pseudo_family.value
                try:
                    self.code_group.selected_code = builder.scf.pw['code'].full_label
                except traitlets.TraitError:
                    # It's not considered a critical issue if the exact code cannot
                    # be retrieved, however we set it to "None" in the interface to
                    # indicate that to the user.
                    self.code_group.selected_code = None
                self.input_structure = process_node.inputs.structure

    def _refresh_outputs_keys(self, process_id):
        with self.hold_trait_notifications():
            process_node = super()._refresh_outputs_keys(process_id)
            if process_node is not None:
                if 'band_structure' in process_node.outputs:
                    self.band_structure = process_node.outputs.band_structure
                elif process_node.is_sealed:
                    self.state = WizardApp.State.FAIL

    def submit(self, _=None):
        assert self.input_structure is not None

        builder = WorkflowFactory('quantumespresso.pw.bands').get_builder()

        builder.scf.pw.code = self.code_group.selected_code
        builder.scf.pw.parameters = load_default_parameters()
        builder.scf.pw.metadata.options = self.options
        builder.scf.kpoints_distance = Float(0.8)
        builder.scf.pseudo_family = Str(self.pseudo_family_selection.value)

        builder.bands.pw.code = self.code_group.selected_code
        builder.bands.pw.parameters = load_default_parameters()
        builder.bands.pw.metadata.options = self.options
        builder.bands.pseudo_family = Str(self.pseudo_family_selection.value)

        builder.structure = self.input_structure

        self.process = submit(builder)

    def reset(self):
        self.process = None
