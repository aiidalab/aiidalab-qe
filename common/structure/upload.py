from __future__ import print_function

from __future__ import absolute_import
import ase.io
import ipywidgets as ipw
from fileupload import FileUploadWidget
import tempfile
import nglview
from six.moves import zip

def get_example_structure(key):
    from ase.io import read
    return read('structures/' + key)


class StructureUploadWidget(ipw.VBox):

    DATA_FORMATS = ('StructureData', 'CifData')

    def __init__(self, text="Upload Structure", **kwargs):
        """ Upload a structure and store it in AiiDA database.

        :param text: Text to display before upload button
        :type text: str
        """

        self.file_upload = FileUploadWidget(text)
        structures = {
                "Select structure": False,
                }
        self.structure_select = ipw.Dropdown(
                options=[],
                description='Or choose from examples:',
                style={'description_width': '160px'},
                disabled=False)
        self.viewer = nglview.NGLWidget()
        self.btn_store = ipw.Button(
            description='Store in AiiDA', disabled=True)
        self.structure_description = ipw.Text(
            placeholder="Description (optional)")

        self.structure_ase = None
        select = ipw.HBox([self.file_upload, self.structure_select])
        store = ipw.HBox([self.btn_store, self.structure_description])
        children = [select, self.viewer, store]

        super(StructureUploadWidget, self).__init__(
            children=children, **kwargs)

        self.file_upload.observe(self._on_file_upload, names='data')
        self.structure_select.observe(self._on_structure_select, names=['value'])
        self.btn_store.on_click(self._on_click_store)

        from aiida import load_dbenv, is_dbenv_loaded
        from aiida.backends import settings
        if not is_dbenv_loaded():
            load_dbenv(profile=settings.AIIDADB_PROFILE)

    # pylint: disable=unused-argument
    def _on_file_upload(self, change):
        self.tmp_folder = tempfile.mkdtemp()
        tmp = self.tmp_folder + '/' + self.file_upload.filename
        with open(tmp, 'w') as f:
            f.write(self.file_upload.data)
        structure_ase = self.get_ase(self.tmp_folder + '/' + self.file_upload.filename)
        self.select_structure(s=structure_ase, name=self.file_upload.filename)

    def _on_structure_select(self, change):
        global atoms
        indx = change['owner'].index
        atoms = change['new']
        if atoms is False:
            self.select_structure(s=None, name=None)
            return None
        formula = atoms.get_chemical_formula()
        self.select_structure(s=atoms, name=formula)


    def select_structure(self, s, name):
        self.btn_store.disabled = False
        if s is None:
            self.structure_ase = None
            self.btn_store.disabled = True
            self.structure_description.value = ""
            self.refresh_view()
            return

        self.structure_description.value = self.get_description(
            s, name)
        self.structure_ase = s
        self.refresh_view()

    def get_ase(self, fname):
        try:
            traj = ase.io.read(fname, index=":")
        except AttributeError:
            print("Looks like {} file does not contain structure coordinates".
                  format(fname))
            return None
        if len(traj) > 1:
            print(
                "Warning: Uploaded file {} contained more than one structure. I take the first one."
                .format(fname))
        return traj[0]

    def get_description(self, structure_ase, name):
        formula = structure_ase.get_chemical_formula()
        return "{} ({})".format(formula, name)

    def refresh_view(self):
        viewer = self.viewer
        # Note: viewer.clear() only removes the 1st component
        # pylint: disable=protected-access
        for comp_id in viewer._ngl_component_ids:
            viewer.remove_component(comp_id)

        if self.structure_ase is None:
            return

        viewer.add_component(nglview.ASEStructure(
            self.structure_ase))  # adds ball+stick
        viewer.add_unitcell()

    # pylint: disable=unused-argument
    def _on_click_store(self, change):
        self.store_structure(
            self.file_upload.filename,
            description=self.structure_description.value)

    def store_structure(self, name, description=None):
        structure_ase = self.structure_ase
        if structure_ase is None:
            return

        from aiida.orm.data.structure import StructureData
        self.structure_node = StructureData(ase=structure_ase)
        if description is None:
            self.structure_node.description = self.get_description(
                structure_ase, name)
        else:
            self.structure_node.description = description
        self.structure_node.label = ".".join(name.split('.')[:-1])
        self.structure_node.store()
        print("Stored in AiiDA: " + repr(self.structure_node))

#EOF
