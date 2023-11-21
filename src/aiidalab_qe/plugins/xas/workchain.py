import os
import tarfile

import requests
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

XspectraCrystalWorkChain = WorkflowFactory("quantumespresso.xspectra.crystal")


def load_or_create_group(group_label_string):
    """Check for the existence of a group matching a label, and create one if not found."""

    from aiida.orm import UpfFamily

    try:
        group = orm.load_group(group_label_string)
    except BaseException:
        group = UpfFamily(group_label_string)
        group.store()

    return group


def check_ch_pseudos_for_elements(group_label):
    """Check a set of core-hole pseudos for which elements are available."""

    from aiida import orm

    elements_list = []
    group = orm.load_group()
    for node in group.nodes:
        element = node.label.split(".")[0]
        elements_list.append(element)

    return elements_list


base_url = "https://github.com/PNOGillespie/Core_Level_Spectra_Pseudos/raw/main"
head_path = "/home/jovyan/Utils/QE/Pseudos"
dir_header = "cls_pseudos"
functionals = ["pbe"]
core_wfc_dir = "core_wfc_data"
gipaw_dir = "gipaw_pseudos"
ch_pseudo_dir = "ch_pseudos/star1s"

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
        remote_archive_filename = f"{base_url}/{func}/{archive_filename}"
        local_archive_filename = f"{dir}/{archive_filename}"

        env = os.environ.copy()
        env["PATH"] = f"{env['PATH']}:{dir}"

        response = requests.get(remote_archive_filename, timeout=30)
        response.raise_for_status()
        with open(local_archive_filename, "wb") as handle:
            handle.write(response.content)
            handle.flush()

        with tarfile.open(local_archive_filename, "r:gz") as tarfil:
            tarfil.extractall(dir)

cls_group_head = load_or_create_group(dir_header)
for func in functionals:
    func_group = load_or_create_group(f"{dir_header}/{func}")
    # Strictly speaking, this one won't have UPF data in it, but it *is* related
    # to the other UPF data anyway.
    core_wfc_group = load_or_create_group(f"{dir_header}/{func}/{core_wfc_dir}")
    core_wfc_node_labels = [node.label for node in core_wfc_group.nodes]
    core_wfc_nodes = []
    core_wfc_path = f"{head_path}/{core_wfc_group.label}"
    for file in os.listdir(core_wfc_path):
        if file not in core_wfc_node_labels:
            new_singlefile = orm.SinglefileData(
                f"{core_wfc_path}/{file}", filename="stdout"
            )
            new_singlefile.label = file
            new_singlefile.store()
            core_wfc_nodes.append(new_singlefile)
    if len(core_wfc_nodes) > 0:
        core_wfc_group.add_nodes(core_wfc_nodes)

    gipaw_group = load_or_create_group(f"{dir_header}/{func}/{gipaw_dir}")
    gipaw_node_labels = [node.label for node in gipaw_group.nodes]
    gipaw_nodes = []
    gipaw_path = f"{head_path}/{gipaw_group.label}"
    for file in os.listdir(gipaw_path):
        if file not in gipaw_node_labels:
            new_upf = orm.UpfData(f"{gipaw_path}/{file}", filename=file)
            new_upf.label = file
            new_upf.store()
            gipaw_nodes.append(new_upf)
    if len(gipaw_nodes) > 0:
        gipaw_group.add_nodes(gipaw_nodes)

    ch_group = load_or_create_group(f"{dir_header}/{func}/{ch_pseudo_dir}")
    ch_node_labels = [node.label for node in ch_group.nodes]
    ch_nodes = []
    ch_path = f"{head_path}/{ch_group.label}"
    for file in os.listdir(ch_path):
        if file not in ch_node_labels:
            new_upf = orm.UpfData(f"{ch_path}/{file}", filename=file)
            new_upf.label = file
            new_upf.store()
            ch_nodes.append(new_upf)
    if len(ch_nodes) > 0:
        ch_group.add_nodes(ch_nodes)


def get_builder(codes, structure, parameters, **kwargs):
    from copy import deepcopy

    protocol = parameters["workchain"]["protocol"]
    xas_parameters = parameters["xas"]
    gipaw_pseudo_group = orm.load_group(xas_parameters["gipaw_pseudo_group"])
    ch_pseudo_group = orm.load_group(xas_parameters["ch_pseudo_group"])
    core_wfc_data_group = orm.load_group(xas_parameters["core_wfc_data_group"])
    # set pseudo for element
    pseudos = {}
    core_wfc_data = {}
    core_hole_treatments = xas_parameters["core_hole_treatments"]
    elements_list = xas_parameters["elements_list"]
    for element in elements_list:
        gipaw_pseudo = [
            gipaw_upf
            for gipaw_upf in gipaw_pseudo_group.nodes
            if gipaw_upf.element == element
        ][0]
        ch_pseudo = [
            ch_upf for ch_upf in ch_pseudo_group.nodes if ch_upf.element == element
        ][0]
        core_wfc_node = [
            sf for sf in core_wfc_data_group.nodes if sf.label.split(".")[0] == element
        ][0]
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
            },
        }
    }
    # TODO: We need to ensure that the family used for selecting the pseudopotentials in
    # the CrystalWorkChain is set to the same as the functional needed for the
    # core-hole/gipaw pseudo pair. In the future, this should check which one it is
    # and set the correct value if it isn't already.
    # For now, I've tried to override the pseudo family automatically using the example below,
    # but it seems that the override is being *overriden* by some other part of the App
    # chosen_pseudo_family = overrides["core"]["scf"]["pseudo_family"]
    # if "/PBE/" not in chosen_pseudo_family:
    #     if chosen_pseudo_family.split("/")[0] == "PseudoDojo":
    #         overrides["core"]["scf"]["pseudo_family"] = "PseudoDojo/0.4/PBE/SR/standard/upf"
    #     elif chosen_pseudo_family.split("/")[0] == "SSSP":
    #         overrides["core"]["scf"]["pseudo_family"] = "SSSP/1.2/PBE/efficiency"

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
