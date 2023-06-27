def test_structure_aiida_database():
    """Test import structure from example."""
    from aiida.plugins import DataFactory
    from ase.io import read

    from aiidalab_qe.app.structures import (
        Examples,
        StructureSelectionStep,
        structure_manager_widget,
    )

    wg = StructureSelectionStep(manager=structure_manager_widget)
    # select structure
    StructureData = DataFactory("core.structure")
    mol = read(Examples[0][1])
    mol = StructureData(ase=mol)
    mol.store()
    #
    structure = wg.manager.children[0].children[2]
    structure.search()
    key = [key for key in structure.results.options if key.startswith(f"PK: {mol.pk}")][
        0
    ]
    structure.results.value = structure.results.options[key]
    wg.confirm()
    assert wg.structure.get_formula() == "Si2"


def test_structure_from_example():
    """Test import structure from example."""

    from aiidalab_qe.app.structures import (
        StructureSelectionStep,
        structure_manager_widget,
    )

    wg = StructureSelectionStep(manager=structure_manager_widget)
    # select structure from example
    structure = wg.manager.children[0].children[3]
    structure.children[0].value = structure.children[0].options[1][1]
    wg.confirm()
    assert wg.structure.get_formula() == "Si2"
