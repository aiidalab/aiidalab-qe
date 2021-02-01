import traitlets

from aiida.orm import StructureData, Float, Str
from aiida.plugins import WorkflowFactory
from aiida.engine import submit

from wizard import WizardApp
from codes import CodeSubmitWidget
from util import load_default_parameters


class RelaxSubmitWidget(CodeSubmitWidget):

    input_structure = traitlets.Instance(StructureData, allow_none=True)
    output_structure = traitlets.Instance(StructureData, allow_none=True)

    def _update_state(self):
        if self.process is None:
            if self.input_structure is None:
                if self.output_structure is None:
                    self.state = WizardApp.State.INIT
                else:
                    self.state = WizardApp.State.FAIL
            else:
                if self.output_structure is None:
                    self.state = WizardApp.State.READY
                else:
                    self.state = WizardApp.State.SUCCESS
        else:
            super()._update_state()

    @traitlets.observe("input_structure")
    def _observe_input_structure(self, change):
        self._update_state()

    @traitlets.observe("output_structure")
    def _observe_output_structure(self, change):
        self._update_state()

    @traitlets.observe("process")
    def _observe_process(self, change):
        with self.hold_trait_notifications():
            super()._observe_process(change)
            process_node = change["new"]
            process_node = super()._observe_process(change)
            if process_node is not None:
                self.input_structure = process_node.inputs.structure
        return process_node

    def skip(self, _):
        self.output_structure = self.input_structure

    def _refresh_outputs_keys(self, process_id):
        process_node = super()._refresh_outputs_keys(process_id)
        if process_node is not None:
            with self.hold_trait_notifications():
                if "output_structure" in process_node.outputs:
                    self.output_structure = process_node.outputs.output_structure
                elif process_node.is_sealed:
                    self.state = WizardApp.State.FAIL

    def submit(self, _=None):
        assert self.input_structure is not None

        builder = WorkflowFactory("quantumespresso.pw.relax").get_builder()
        builder.base.pw.code = self.code_group.selected_code
        builder.base.pw.parameters = load_default_parameters()
        builder.base.pw.metadata.options = self.options
        builder.base.kpoints_distance = Float(0.8)
        builder.base.pseudo_family = Str(self.pseudo_family_selector.value)
        builder.structure = self.input_structure

        self.process = submit(builder)

    def reset(self):
        with self.hold_trait_notifications():
            self.process = None
            self.output_structure = None
