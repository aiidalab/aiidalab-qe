import os
import tarfile
from importlib import resources
from pathlib import Path

import requests
import yaml
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

from aiidalab_qe.plugins import xas as xas_folder

XspectraCrystalWorkChain = WorkflowFactory("quantumespresso.xspectra.crystal")
PSEUDO_TOC = yaml.safe_load(resources.read_text(xas_folder, "pseudo_toc.yaml"))
pseudo_data_dict = PSEUDO_TOC["pseudos"]
xch_elements = PSEUDO_TOC["xas_xch_elements"]

base_url = "https://github.com/PNOGillespie/Core_Level_Spectra_Pseudos/raw/main"
head_path = f"{Path.home()}/.local/lib"
dir_header = "cls_pseudos"
functionals = ["pbe"]
core_wfc_dir = "core_wfc_data"
gipaw_dir = "gipaw_pseudos"
ch_pseudo_dir = "ch_pseudos/star1s"


def _load_or_import_nodes_from_filenames(in_dict, path, core_wfc_data=False):
    for filename in in_dict.values():
        try:
            orm.load_node(filename)
        except BaseException:
            if not core_wfc_data:
                new_upf = orm.UpfData(f"{path}/{filename}", filename=filename)
                new_upf.label = filename
                new_upf.store()
            else:
                new_singlefile = orm.SinglefileData(
                    f"{path}/{filename}", filename="stdout"
                )
                new_singlefile.label = filename
                new_singlefile.store()


def _download_extract_pseudo_archive(func):
    dir = f"{head_path}/{dir_header}/{func}"
    archive_filename = f"{func}_ch_pseudos.tgz"
    remote_archive_filename = f"{base_url}/{func}/{archive_filename}"
    local_archive_filename = f"{dir}/{archive_filename}"

    env = os.environ.copy()
    env["PATH"] = f"{env['PATH']}:{Path.home() / '.local' / 'lib'}"

    response = requests.get(remote_archive_filename, timeout=30)
    response.raise_for_status()
    with open(local_archive_filename, "wb") as handle:
        handle.write(response.content)
        handle.flush()
        response.close()

    with tarfile.open(local_archive_filename, "r:gz") as tarfil:
        tarfil.extractall(dir)


url = f"{base_url}"
for func in functionals:
    dir = f"{head_path}/{dir_header}/{func}"
    os.makedirs(dir, exist_ok=True)
    archive_filename = f"{func}_ch_pseudos.tgz"
    archive_found = False
    for entry in os.listdir(dir):
        if entry == archive_filename:
            archive_found = True
    if not archive_found:
        _download_extract_pseudo_archive(func)


# Check all the pseudos/core-wfc data files in the TOC dictionary
# and load/check all of them before proceeding. Note that this
# approach relies on there not being multiple instances of nodes
# with the same label.
for func in functionals:
    gipaw_pseudo_dict = pseudo_data_dict[func]["gipaw_pseudos"]
    core_wfc_dict = pseudo_data_dict[func]["core_wavefunction_data"]
    core_hole_pseudo_dict = pseudo_data_dict[func]["core_hole_pseudos"]
    main_path = f"{head_path}/{dir_header}/{func}"
    core_wfc_dir = f"{main_path}/core_wfc_data"
    gipaw_dir = f"{main_path}/gipaw_pseudos"
    ch_pseudo_dir = f"{main_path}/ch_pseudos/star1s"
    # First, check that the local directories contain what's in the pseudo_toc
    for pseudo_dir, pseudo_dict in zip(
        [gipaw_dir, core_wfc_dir, ch_pseudo_dir],
        [gipaw_pseudo_dict, core_wfc_dict, core_hole_pseudo_dict],
    ):
        pseudo_toc_mismatch = os.listdir(pseudo_dir) != pseudo_dict.values()

    # Re-download the relevant archive if there is a mismatch
    if pseudo_toc_mismatch:
        _download_extract_pseudo_archive(func)

    _load_or_import_nodes_from_filenames(
        in_dict=gipaw_pseudo_dict,
        path=gipaw_dir,
    )
    _load_or_import_nodes_from_filenames(
        in_dict=core_wfc_dict, path=core_wfc_dir, core_wfc_data=True
    )
    _load_or_import_nodes_from_filenames(
        in_dict=core_hole_pseudo_dict["1s"], path=ch_pseudo_dir
    )


def get_builder(codes, structure, parameters, **kwargs):
    from copy import deepcopy

    protocol = parameters["workchain"]["protocol"]
    xas_parameters = parameters["xas"]
    gipaw_pseudos_dict = pseudo_data_dict["pbe"]["gipaw_pseudos"]
    ch_pseudos_dict = pseudo_data_dict["pbe"]["core_hole_pseudos"]["1s"]
    core_wfc_data_dict = pseudo_data_dict["pbe"]["core_wavefunction_data"]
    # set pseudo for element
    pseudos = {}
    core_wfc_data = {}
    core_hole_treatments = xas_parameters["core_hole_treatments"]
    elements_list = xas_parameters["elements_list"]
    for element in elements_list:
        gipaw_pseudo = orm.load_node(gipaw_pseudos_dict[element])
        ch_pseudo = orm.load_node(ch_pseudos_dict[element])
        core_wfc_node = orm.load_node(core_wfc_data_dict[element])
        pseudos[element] = {"gipaw": gipaw_pseudo, "core_hole": ch_pseudo}
        core_wfc_data[element] = core_wfc_node

    # TODO should we override the cutoff_wfc, cutoff_rho by the new pseudo?
    # In principle we should, if we know what that value is, but that would
    # require testing them first...

    # (13/10/23) I'm keeping the part about molecules in for future reference,
    # but we need to establish the protocol & backend code for XAS of molecules
    # before thinking about a workflow.
    is_molecule_input = (
        True if xas_parameters.get("structure_type") == "molecule" else False
    )

    # core_hole_treatment = xas_parameters["core_hole_treatment"]
    # core_hole_treatments = {element: core_hole_treatment for element in elements_list}

    structure_preparation_settings = {
        # "supercell_min_parameter": Float(supercell_min_parameter_map[protocol]),
        "is_molecule_input": orm.Bool(is_molecule_input),
    }
    spglib_settings = orm.Dict({"symprec": 1.0e-3})

    pw_code = codes["pw"]
    xs_code = codes["xspectra"]
    overrides = {
        "core": {
            "scf": deepcopy(parameters["advanced"]),
            # PG: Here, we set a "variable" broadening scheme, which actually defines a constant broadening
            # The reason for this is that in "gamma_mode = constant", the Lorenzian broadening parameter
            # is defined by "xgamma" (in "PLOT"), but this parameter *also* controls the broadening value
            # used in the Lanczos algorithm to enhance the convergence rate. In order to give the user a
            # final spectrum with minimal broadening, we use "gamma_mode = variable", which uses a different
            # parameter set ("gamma_energy(1-2)", "gamma_value(1-2)") and thus allows us to decouple spectrum
            # broadening from Lanczos broadening and avoid having to re-plot the final spectrum.
            "xs_prod": {
                "xspectra": {
                    "parameters": {
                        "PLOT": {
                            "gamma_mode": "variable",
                            "gamma_energy(1)": 0,
                            "gamma_energy(2)": 1,
                            "gamma_value(1)": 0.1,
                            "gamma_value(2)": 0.1,
                        }
                    }
                }
            }
        }
    }

    builder = XspectraCrystalWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        xs_code=xs_code,
        structure=structure,
        protocol=protocol,
        pseudos=pseudos,
        elements_list=elements_list,
        core_hole_treatments=core_hole_treatments,
        core_wfc_data=core_wfc_data,
        structure_preparation_settings=structure_preparation_settings,
        electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
        spin_type=SpinType(parameters["workchain"]["spin_type"]),
        # TODO: We will need to merge the changes in AiiDA-QE PR#969 in order
        # to better handle magnetic and Hubbard data. For now, we can probably
        # leave it as it is.
        initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
        overrides=overrides,
        **kwargs,
    )
    builder.pop("relax")
    builder.pop("clean_workdir", None)
    builder.spglib_settings = spglib_settings
    # there is a bug in aiida-quantumespresso xps, that one can not set the kpoints
    # this is fxied in a PR, but we need to wait for the next release.
    # we set a large kpoints_distance value to set the kpoints to 1x1x1
    if is_molecule_input:
        # kpoints = KpointsData()
        # kpoints.set_kpoints_mesh([1, 1, 1])
        # parameters["advanced"]["kpoints"] = kpoints
        # builder.ch_scf.kpoints_distance = Float(5)
        pass
    return builder


workchain_and_builder = {
    "workchain": XspectraCrystalWorkChain,
    "exclude": ("clean_workdir", "structure", "relax"),
    "get_builder": get_builder,
}
