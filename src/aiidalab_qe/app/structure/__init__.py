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
)
from aiidalab_qe.common.infobox import InAppGuide
from aiidalab_qe.common.widgets import CategorizedStructureExamplesWidget, QeWizardStep
from aiidalab_widgets_base import (
    BasicCellEditor,
    BasicStructureEditor,
    StructureManagerWidget,
    StructureUploadWidget,
)

# The Examples list of (name, file) tuple curretly passed to
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


class StructureSelectionStep(QeWizardStep[StructureStepModel]):
    """Integrated widget for the selection and edition of structure.
    The widget includes a structure manager that allows to select a structure
    from different sources. It also includes the structure editor. Both the
    structure importers and the structure editors can be extended by plugins.
    """

    def __init__(self, model: StructureStepModel, **kwargs):
        super().__init__(model=model, **kwargs)
        self._model.observe(
            self._on_confirmation_change,
            "confirmed",
        )
        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )

    def _render(self):
        examples_by_category = {"Simple crystals": Examples}
        plugin_structure_examples = {
            item["title"]: item["structures"]
            for item in get_entry_items(
                "aiidalab_qe.properties", "structure_examples"
            ).values()
        }
        examples_by_category.update(plugin_structure_examples)

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
            AddingTagsEditor(title="Edit StructureData"),
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

        self.confirm_button = ipw.Button(
            description="Confirm",
            tooltip="Confirm the currently selected structure and go to the next step.",
            button_style="success",
            icon="check-circle",
            layout=ipw.Layout(width="auto"),
        )
        ipw.dlink(
            (self, "state"),
            (self.confirm_button, "disabled"),
            lambda state: state != self.State.CONFIGURED,
        )
        self.confirm_button.on_click(self.confirm)

        self.message_area = ipw.HTML()
        ipw.dlink(
            (self._model, "message_area"),
            (self.message_area, "value"),
        )

        self.children = [
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
            self.message_area,
            self.confirm_button,
        ]

    def is_saved(self):
        return self._model.confirmed

    def confirm(self, _=None):
        self.manager.store_structure()
        self._model.message_area = ""
        self._model.confirm()

    def can_reset(self):
        return self._model.confirmed

    def reset(self):
        self._model.reset()

    def _on_input_structure_change(self, _):
        self._model.update_widget_text()
        self._update_state()

    def _on_confirmation_change(self, _):
        self._update_state()

    def _update_state(self):
        if self._model.confirmed:
            self.state = self.State.SUCCESS
        elif self._model.input_structure is None:
            self.state = self.State.READY
        else:
            self.state = self.State.CONFIGURED
