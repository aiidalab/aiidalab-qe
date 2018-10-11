from numpy import array
from aiida.orm import DataFactory
from aiida.work.workfunction import workfunction

StructureData = DataFactory('structure')

@workfunction
def create_diamond_fcc(element, alat):
    """
    Workfunction to create a diamond crystal structure with a given element.

    :param element: The element to create the structure with.
    :param alat: The lattice parameter in Angstrom
    :return: The constructed StructureData object
    """
    the_cell = array([[0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]]) * alat
    structure = StructureData(cell=the_cell)
    structure.append_atom(position=(0., 0., 0.), symbols=str(element))
    structure.append_atom(position=(0.25 * alat, 0.25 * alat, 0.25 * alat), symbols=str(element))

    return structure


@workfunction
def scale_structure(structure, scaling_factor):
    """
    Workfunction to scale a structure

    :param structure: An AiiDA structure to scale
    :param scaling_factor: The scaling factor
    :return: The scaled StructureData object
    """
    ase_old = structure.get_ase()
    ase_new = ase_old.copy()
    ase_new.set_cell(ase_old.get_cell() * float(scaling_factor), scale_atoms=True)
    structure = StructureData(ase=ase_new)

    return structure