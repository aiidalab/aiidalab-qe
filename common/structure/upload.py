from __future__ import print_function

from __future__ import absolute_import
import ase.io
import ipywidgets as ipw
from fileupload import FileUploadWidget
import tempfile
import nglview
from six.moves import zip

def get_example_structure(key):
    from ase.lattice.spacegroup import crystal
    from ase import Atoms

    if key == 'diamond':
        # This is the lattice constant in angstrom
        alat = 3.56
        diamond_ase = crystal('C', [(0,0,0)], spacegroup=227,
                          cellpar=[alat, alat, alat, 90, 90, 90],primitive_cell=True)
        return diamond_ase
    elif key == 'al':
        # This is the lattice constant in angstrom
        alat = 4.05
        Al_ase = crystal('Al', [(0,0,0)], spacegroup=225,
                          cellpar=[alat, alat, alat, 90, 90, 90],primitive_cell=False)
        return Al_ase
    elif key == 'si':
        cell = [[2.6954645, 2.6954645, 0],
                [2.6954645, 0, 2.6954645],
                [0, 2.6954645, 2.6954645]]
        a_si = Atoms("Si2", cell=cell, scaled_positions=[[0, 0, 0], [0.25, 0.25, 0.25]])

        alat = 2.69
        a_si = crystal('Si', [(0, 0, 0)], spacegroup=227,
                cellpar=[alat, alat, alat, 90, 90, 90],primitive_cell=True)
        return a_si
    elif key == 'gaas':
        # This is the lattice constant in angstrom
        alat = 5.75
        GaAs_ase = crystal('GaAs', [(0,0,0),(0.25,0.25,0.25)], spacegroup=216,
                          cellpar=[alat, alat, alat, 90, 90, 90],primitive_cell=True)
        return GaAs_ase
    elif key == 'co':
        # These are the lattice constants in angstrom
        a = 2.5
        c = 4.07
        Co_ase = crystal('Co', [(1./3,2./3,0.25)], spacegroup=194,
                          cellpar=[a, a, c, 90, 90, 120],primitive_cell=True)
        return Co_ase
    else:
        raise ValueError("Unknown or unsupported example structure '{}'".format(key))




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
#                'Aluminium':get_example_structure('al'),
                'Silicon' : get_example_structure('si'),
                'Diamond':get_example_structure('diamond'),
                'Gallium arsenide': get_example_structure('gaas'),
                }
        self.structure_select = ipw.Dropdown(
                options=structures,
                value=False,
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
        self.structure_select.observe(self._on_structure_select, names='value')
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
        structure_ase = self.get_ase(self.tmp_folder + '/' + name)
        if structure_ase is None:
            return

        from aiida.orm.data.structure import StructureData
        structure_node = StructureData(ase=structure_ase)
        if description is None:
            structure_node.description = self.get_description(
                structure_ase, name)
        else:
            structure_node.description = description
        structure_node.label = ".".join(name.split('.')[:-1])
        structure_node.store()
        print("Stored in AiiDA: " + repr(structure_node))

