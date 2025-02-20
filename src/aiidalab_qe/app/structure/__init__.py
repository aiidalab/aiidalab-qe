"""Widgets for the upload and selection of structure data.

Authors: AiiDAlab team
"""

import pathlib

import ipywidgets as ipw

from aiidalab_qe.app.structure.model import StructureStepModel
from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common import (
    AddingTagsEditor,
    LazyLoadedOptimade,
    LazyLoadedStructureBrowser,
    PeriodicityEditor,
    ShakeNBreakEditor,
)
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.widgets import CategorizedStructureExamplesWidget
from aiidalab_qe.common.wizard import QeConfirmableWizardStep
from aiidalab_widgets_base import (
    BasicCellEditor,
    BasicStructureEditor,
    StructureManagerWidget,
    StructureUploadWidget,
)

# The Examples list of (name, file) tuple currently passed to
# StructureExamplesWidget.
file_path = pathlib.Path(__file__).parent
Examples = [
    ("Bulk silicon", file_path / "examples" / "Si.cif"),
    ("Diamond", file_path / "examples" / "Diamond.cif"),
    ("Gallium arsenide", file_path / "examples" / "GaAs.cif"),
    ("Gold", file_path / "examples" / "Au.cif"),
    ("Nickel", file_path / "examples" / "Ni.cif"),
    ("LiCoO2", file_path / "examples" / "LiCoO2.cif"),
    ("Silicon oxide (alpha quartz)", file_path / "examples" / "SiO2.cif"),
]


class StructureSelectionStep(QeConfirmableWizardStep[StructureStepModel]):
    """Integrated widget for the selection and edition of structure.
    The widget includes a structure manager that allows to select a structure
    from different sources. It also includes the structure editor. Both the
    structure importers and the structure editors can be extended by plugins.
    """

    def __init__(self, model: StructureStepModel, auto_setup=True, **kwargs):
        super().__init__(
            model=model,
            confirm_kwargs={
                "tooltip": "Confirm the currently selected structure and go to the next step",
            },
            **kwargs,
        )
        self._model.observe(
            self._on_installation_change,
            ["installing_sssp", "sssp_installed"],
        )
        self._model.observe(
            self._on_sssp_installed,
            "sssp_installed",
        )
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._install_sssp(auto_setup)

    def _render(self):
        super()._render()

        examples_by_category = {"Simple crystals": Examples}
        plugin_structure_examples = {
            item["title"]: item["structures"]
            for item in get_entry_items(
                "aiidalab_qe.properties", "structure_examples"
            ).values()
        }
        examples_by_category |= plugin_structure_examples

        importers = [
            StructureUploadWidget(title="Upload file"),
            LazyLoadedOptimade(title="OPTIMADE"),
            LazyLoadedStructureBrowser(title="AiiDA database"),
            CategorizedStructureExamplesWidget(
                title="From examples", examples_by_category=examples_by_category
            ),
        ]

        plugin_importers = get_entry_items("aiidalab_qe.properties", "importer")
        importers.extend([importer() for importer in plugin_importers.values()])

        editors = [
            BasicCellEditor(title="Edit cell"),
            BasicStructureEditor(title="Edit structure"),
            AddingTagsEditor(title="Edit atom tags"),
            PeriodicityEditor(title="Edit periodicity"),
            ShakeNBreakEditor(title="ShakeNBreak"),
        ]

        plugin_editors = get_entry_items("aiidalab_qe.properties", "editor")
        editors.extend([editor() for editor in plugin_editors.values()])

        self.manager = StructureManagerWidget(
            importers=importers,
            editors=editors,
            node_class="StructureData",
            storable=False,
            configuration_tabs=[
                "Cell",
                "Selection",
                "Appearance",
                "Download",
            ],
        )

        if self._model.confirmed:  # loaded from a process
            # NOTE important to do this prior to setting up the links
            # to avoid an override of the structure in the model,
            # which in turn would trigger a reset of the model
            self.manager.input_structure = self._model.input_structure

        ipw.dlink(
            (self.manager, "structure_node"),
            (self._model, "input_structure"),
        )
        ipw.link(
            (self._model, "manager_output"),
            (self.manager.output, "value"),
        )

        self.structure_name_text = ipw.Text(
            placeholder="[No structure selected]",
            description="Selected:",
            disabled=True,
            layout=ipw.Layout(width="auto", flex="1 1 auto"),
        )
        ipw.dlink(
            (self._model, "structure_name"),
            (self.structure_name_text, "value"),
        )

        self.content.children = [
            InAppGuide(identifier="structure-step"),
            ipw.HTML("""
                <p>
                    Select a structure from one of the following sources, then
                    click <span style="color: #4caf50;">
                        <i class="fa fa-check-circle"></i> <b>Confirm</b>
                    </span> to go to the next step
                </p>
            """),
            self.manager,
            self.structure_name_text,
        ]

        self.confirm_box.children += (self.sssp_installation,)

        self.children = [
            self.content,
            self.confirm_box,
        ]

    def _post_render(self):
        # After rendering the widget, nglview needs to be resized
        # to properly display the structure
        self.manager.viewer._viewer.handle_resize()

    def confirm(self, _=None):
        self.manager.store_structure()
        super().confirm()

    def can_reset(self):
        return self._model.confirmed

    def reset(self):
        self._model.reset()

    def _on_installation_change(self, _):
        self._model.update_blockers()

    def _on_sssp_installed(self, _):
        self._toggle_sssp_installation_widget()

    def _on_input_structure_change(self, _):
        self._model.update_widget_text()
        self._update_state()

    def _install_sssp(self, auto_setup):
        from aiidalab_qe.common.setup_pseudos import PseudosInstallWidget

        self.sssp_installation = PseudosInstallWidget(auto_start=False)
        ipw.dlink(
            (self.sssp_installation, "busy"),
            (self._model, "installing_sssp"),
        )
        ipw.dlink(
            (self.sssp_installation, "installed"),
            (self._model, "installing_sssp"),
            lambda installed: not installed,
        )
        ipw.dlink(
            (self.sssp_installation, "installed"),
            (self._model, "sssp_installed"),
        )
        if auto_setup:
            self.sssp_installation.refresh()

    def _toggle_sssp_installation_widget(self):
        sssp_installation_display = "none" if self._model.sssp_installed else "block"
        self.sssp_installation.layout.display = sssp_installation_display

    def _update_state(self):
        if self._model.confirmed:
            self.state = self.State.SUCCESS
        elif self._model.input_structure is None:
            self.state = self.State.READY
        else:
            self.state = self.State.CONFIGURED
