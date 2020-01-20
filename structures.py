"""Widgets for the upload and selection of structure data.

Authors:

    * Carl Simon Adorf <simon.adorf@epfl.ch>
"""
import os
import sys
import tempfile

import ase
import nglview
import traitlets
import ipywidgets as ipw

from wizard import WizardApp, WizardAppStep


def _read_ase_atoms(data, filename=None):
    """Read the structure data into an ase.Atoms object."""
    ext = None if filename is None else os.path.splitext(filename)[1]

    try:
        with tempfile.NamedTemporaryFile(suffix=ext) as file:
            file.write(data)
            file.flush()

            traj = ase.io.read(file.name, index=':')
        # TODO: REACT TO len(traj) > 1
    except Exception as error:
        print(error, file=sys.stderr)
    else:
        return traj[0]


class StructureFileUploadWidget(ipw.VBox):
    """Select a structure from a file upload."""

    structure = traitlets.Instance(ase.atoms.Atoms, allow_none=True)

    def __init__(self):
        self.description = ipw.Label("Upload a structure file:")
        self.supported_formats = ipw.HTML(
            """<a href="https://wiki.fysik.dtu.dk/ase/_modules/ase/io/formats.html" target="_blank">
            Supported formats
            </a>""")
        self.file_upload = ipw.FileUpload(
            description="Select file",
            multiple=False,
        )

        self.file_upload.observe(self.on_file_uploaded, names='value')

        super().__init__(children=[self.description, self.file_upload, self.supported_formats])

    def on_file_uploaded(self, change):
        for fn, item in change['new'].items():
            self.structure = _read_ase_atoms(data=item['content'], filename=fn)
            self.file_upload.value.clear()
            break

    def freeze(self):
        self.file_upload.disabled = True

    def reset(self):
        with self.hold_trait_notifications():
            self.file_upload.value.clear()
            self.file_upload.disabled = False


class SelectionStructureUploadWidget(ipw.Dropdown):
    """Select a structure from a given set of of options."""

    structure = traitlets.Instance(ase.atoms.Atoms, allow_none=True)

    def __init__(self, options=None, hint='Select structure', **kwargs):
        if options is None:
            options = []
        else:
            options.insert(0, (hint, None))

        super().__init__(options=options, **kwargs)

    @traitlets.observe('index')
    def _observe_index(self, change):
        index = change['new']
        if index is None:
            self.structure = None
        else:
            filename = self.options[index][1]
            if filename is None:
                self.structure = None
            else:
                traj = ase.io.read(filename, index=':')
                self.structure = traj[0]

    def reset(self):
        with self.hold_trait_notifications():
            self.index = 0
            self.disabled = False

    def freeze(self):
        self.disabled = True


class StructureUploadComboWidget(ipw.VBox, WizardAppStep):
    """Integrated widget for the selection of structures from different sources."""

    structure = traitlets.Instance(ase.atoms.Atoms, allow_none=True)
    confirmed_structure = traitlets.Instance(ase.atoms.Atoms, allow_none=True)

    def __init__(self, data_importers=None, examples=None, viewer=True, **kwargs):
        if data_importers is None:
            self.data_importers = [('Upload', StructureFileUploadWidget())]
        else:
            self.data_importers = data_importers

        if examples:
            self.example_widget = SelectionStructureUploadWidget(options=examples)
            self.data_importers.append(('Examples', self.example_widget))
        else:
            self.example_widget = None

        if len(self.data_importers) > 1:
            self.structure_sources_tab = ipw.Tab(children=[s[1] for s in self.data_importers])
            for i, source in enumerate(self.data_importers):
                self.structure_sources_tab.set_title(i, source[0])
        else:
            self.structure_sources_tab = self.data_importers[0][1]

        self.structure_sources_tab.layout = ipw.Layout(
            display='flex',
            flex_flow='column',
        )

        if viewer:
            self.viewer = nglview.NGLWidget(width='300px', height='300px')
            self.viewer_box = ipw.Box(children=[self.viewer], layout=ipw.Layout(border='solid 1px'))
        else:
            self.viewer = None
            self.viewer_box = None

        self.structure_name_text = ipw.Text(
            placeholder='[No structure selected]',
            description='Selected:',
            disabled=True,
            layout=ipw.Layout(width='auto', flex="1 1 auto"),
        )

        self.structure_sources_tab.layout = ipw.Layout(min_width='600px')

        grid = ipw.GridspecLayout(1, 3)
        grid[:, :2] = self.structure_sources_tab
        grid[:, 2] = self.viewer_box

        self.confirm_button = ipw.Button(
            description='Confirm',
            button_style='success',
            icon='check-circle',
            disabled=True,
            layout=ipw.Layout(width='auto'),
        )
        self.confirm_button.on_click(self.confirm)

        for data_importer in self.data_importers:
            traitlets.dlink((data_importer[1], 'structure'), (self, 'structure'))

        super().__init__(
            children=[grid, self.structure_name_text, self.confirm_button], **kwargs)

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
                self.structure_name_text.value = str(self.structure.get_chemical_formula())
            self._update_state()
            self.refresh_view()

    @traitlets.observe('confirmed_structure')
    def _observe_confirmed_structure(self, change):
        with self.hold_trait_notifications():
            self._update_state()

    @traitlets.observe('state')
    def _observe_state(self, change):
        with self.hold_trait_notifications():
            state = change['new']
            if state is WizardApp.State.SUCCESS:
                self.freeze()
            self.confirm_button.disabled = self.state != WizardApp.State.CONFIGURED

    def freeze(self):
        for child in self.structure_sources_tab.children:
            child.freeze()

    def refresh_view(self):
        # Note: viewer.clear() only removes the 1st component (TODO: FIX UPSTREAM!)
        for comp_id in self.viewer._ngl_component_ids:
            self.viewer.remove_component(comp_id)
        if self.structure is not None:
            self.viewer.add_component(nglview.ASEStructure(self.structure))
            self.viewer.add_unitcell()

    def confirm(self, button=None):
        with self.hold_trait_notifications():
            self.confirmed_structure = self.structure

    def reset(self):  # unconfirm
        with self.hold_trait_notifications():
            for child in self.structure_sources_tab.children:
                child.reset()
            self.confirmed_structure = None
            self.structure = None
