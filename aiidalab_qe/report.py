from aiida.orm import WorkChainNode
from aiida.plugins import WorkflowFactory

PwBaseWorkChain = WorkflowFactory("quantumespresso.pw.base")


FUNCTIONAL_LINK_MAP = {
    "PBE": "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.77.3865",
    "PBEsol": "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.136406",
}

PSEUDO_LINK_MAP = {
    "SSSP": "https://www.materialscloud.org/discover/sssp/table/efficiency"
}

PROTOCOL_PSEUDO_MAP = {
    "fast": "SSSP/1.2/PBE/efficiency",
    "moderate": "SSSP/1.2/PBE/efficiency",
    "precise": "SSSP/1.2/PBE/precision",
}

FUNCTIONAL_REPORT_MAP = {
    "LDA": "local density approximation (LDA)",
    "PBE": "generalized gradient approximation of Perdew-Burke-Ernzerhof (PBE)",
    "PBEsol": "the revised generalized gradient approximation of Perdew-Burke-Ernzerhof (PBE) for solids",
}

PROTOCOL_BANDS_KPOINTS_DISTANCE = {
    "fast": 0.1,
    "moderate": 0.025,
    "precise": 0.015,
}
def _generate_report_dict(qeapp_wc: WorkChainNode):
    builder_parameters = qeapp_wc.base.extras.get("builder_parameters", {})

    # Properties
    run_relax = builder_parameters.get("relax_type") != "none"
    run_bands = builder_parameters.get("run_bands")
    run_pdos = builder_parameters.get("run_pdos")

    yield "relaxed", run_relax
    yield "relax_method", builder_parameters["relax_type"]
    yield "bands_computed", run_bands
    yield "pdos_computed", run_pdos

    # Material settings
    yield "material_magnetic", builder_parameters["spin_type"]
    yield "electronic_type", builder_parameters["electronic_type"]
    yield "tot_charge", builder_parameters["tot_charge"]
    yield "periodicity", builder_parameters["periodicity"]

    if builder_parameters["periodicity"] == "xy":
        yield "kpoint_path" , builder_parameters["two_dim_kpoints_path"]

    # Calculation settings
    yield "protocol", builder_parameters["protocol"]

    try:
        pseudo_family = builder_parameters.get("pseudo_family", None)
        if pseudo_family is None:
            protocol = builder_parameters.get("pseudo_family", "moderate")
            pseudo_family = PROTOCOL_PSEUDO_MAP[protocol]

        yield "pseudo_family", pseudo_family

        pseudo_family_list = pseudo_family.split("/")
        pseudo_library = pseudo_family_list[0]
        yield "pseudo_library", pseudo_library

        if pseudo_library == "SSSP":
            yield "pseudo_version", pseudo_family_list[1]
            functional = pseudo_family_list[2]
            yield "functional", functional
            yield "pseudo_protocol", pseudo_family_list[3]
        else:
            raise NotImplementedError
    except (KeyError, AttributeError):
        pass
    yield "pseudo_link", PSEUDO_LINK_MAP[pseudo_library]
    yield "functional_link", FUNCTIONAL_LINK_MAP[functional]

    pw_parameters = None
    energy_cutoff_wfc = None
    energy_cutoff_rho = None
    degauss = None
    smearing = None
    scf_kpoints_distance = None
    bands_kpoints_distance = None

    try:
        scf_kpoints_distance = qeapp_wc.inputs.kpoints_distance_override.value
    except AttributeError:
        scf_kpoints_distance = None
    nscf_kpoints_distance = None

    default_params = PwBaseWorkChain.get_protocol_inputs(
        builder_parameters["protocol"]
    )["pw"]["parameters"]["SYSTEM"]
    try:
        degauss = qeapp_wc.inputs.degauss_override.value
    except AttributeError:
        # read default from protocol
        degauss = default_params["degauss"]

    try:
        smearing = qeapp_wc.inputs.smearing_override.value
    except AttributeError:
        # read default from protocol
        smearing = default_params["smearing"]

    # The run_relax variable is not same as if statement below.
    # the "relax" port is poped out to skip the real relaxiation
    # which is not the case of SCF, we use the relax workchain but with
    # relax_type set to none as SCF calculation.
    if "relax" in qeapp_wc.inputs:
        pw_parameters = qeapp_wc.inputs.relax.base.pw.parameters.get_dict()
        if scf_kpoints_distance is None:
            scf_kpoints_distance = qeapp_wc.inputs.relax.base.kpoints_distance.value

    if run_bands:
        pw_parameters = qeapp_wc.inputs.bands.scf.pw.parameters.get_dict()
        if scf_kpoints_distance is None:
            scf_kpoints_distance = qeapp_wc.inputs.bands.scf.kpoints_distance.value
        try: 
            bands_kpoints_distance = qeapp_wc.inputs.bands.bands_kpoints_distance.value
        except AttributeError:
            bands_kpoints_distance = PROTOCOL_BANDS_KPOINTS_DISTANCE[builder_parameters["protocol"]]
    if run_pdos:
        scf_kpoints_distance = (
            scf_kpoints_distance or qeapp_wc.inputs.pdos.scf.kpoints_distance.value
        )
        pw_parameters = (
            pw_parameters or qeapp_wc.inputs.pdos.scf.pw.parameters.get_dict()
        )
        nscf_kpoints_distance = qeapp_wc.inputs.pdos.nscf.kpoints_distance.value

    energy_cutoff_wfc = round(pw_parameters["SYSTEM"]["ecutwfc"])
    energy_cutoff_rho = round(pw_parameters["SYSTEM"]["ecutrho"])

    yield "energy_cutoff_wfc", energy_cutoff_wfc
    yield "energy_cutoff_rho", energy_cutoff_rho
    yield "degauss", degauss
    yield "smearing", smearing
    yield "scf_kpoints_distance", scf_kpoints_distance
    yield "bands_kpoints_distance", bands_kpoints_distance
    yield "nscf_kpoints_distance", nscf_kpoints_distance


def generate_report_dict(qeapp_wc):
    """Generate a dictionary for reporting the inputs for the `QeAppWorkChain`"""
    return dict(_generate_report_dict(qeapp_wc))


def generate_report_text(report_dict):
    """Generate a text for reporting the inputs for the `QeAppWorkChain`

    :param report_dict: dictionary generated by the `generate_report_dict` function.
    """

    report_string = (
        "All calculations are performed within the density-functional "
        "theory formalism as implemented in the Quantum ESPRESSO code. "
        "The pseudopotential for each element is extracted from the "
        f'{report_dict["Pseudopotential library"][0]} '
        "library. The wave functions "
        "of the valence electrons are expanded in a plane wave basis set, using an "
        "energy cutoff equal to "
        f'{round(report_dict["Plane wave energy cutoff (wave functions)"][0])} Ry '
        "for the wave functions and "
        f'{round(report_dict["Plane wave energy cutoff (charge density)"][0])} Ry '
        "for the charge density and potential. "
        "The exchange-correlation energy is "
        "calculated using the "
        f'{FUNCTIONAL_REPORT_MAP[report_dict["Functional"][0]]}. '
        "A Monkhorst-Pack mesh is used for sampling the Brillouin zone, where the "
        "distance between the k-points is set to "
    )
    kpoints_distances = []
    kpoints_calculations = []

    for calc in ("SCF", "NSCF", "Bands"):
        if f"K-point mesh distance ({calc})" in report_dict:
            kpoints_distances.append(
                str(report_dict[f"K-point mesh distance ({calc})"][0])
            )
            kpoints_calculations.append(calc)

    report_string += ", ".join(kpoints_distances)
    report_string += " for the "
    report_string += ", ".join(kpoints_calculations)
    report_string += " calculation"
    if len(kpoints_distances) > 1:
        report_string += "s, respectively"
    report_string += "."

    return report_string


def get_hubbard_occupations_list(qeapp_wc: WorkChainNode) -> list:
    """Get the occupations dictionary for the Hubbard U calculation.

    :param qeapp_wc: `QeAppWorkChain` instance.
    :returns: list of dictionaries of Hubbard occupations.
    """
    builder_parameters = qeapp_wc.base.extras.get("builder_parameters", {}) # To know from where to get the occupations! 
    workchain = qeapp_wc.called[-1]
    
    initial_pattern = 'HUBBARD OCCUPATIONS'
    end_pattern = 'Number of occupied Hubbard levels ='
    list_string = []
    with workchain.called[0].outputs.retrieved.base.repository.open('aiida.out', 'r') as file:
        file_contents = file.readlines()
        i = 0
        while i < len(file_contents):
            if initial_pattern in file_contents[i]:
                i += 1
                text = ''
                while i < len(file_contents) and end_pattern not in file_contents[i]:
                    text += file_contents[i]
                    i += 1
                if text:
                    list_string.append(text)
            else:
                i += 1
          
    last_text = list_string[-1]
    lines = last_text.splitlines()

    eigenvalues_up = []
    eigenvalues_down = []
    occupation_matrix_up = []
    occupation_matrix_down = []
    magnetic_moment = [] # magneton/bohr 
    atom_index = []
    trace_occupation = []
    spin_flag = None

    for index, line in enumerate(lines):
        if 'Tr[ns' in line:
            line_temp = line.split()
            trace_occupation.append([line_temp[-3], line_temp[-2],line_temp[-1]])
        elif 'ATOM' in line:
            atom_index.append(line.split()[-2])
        elif 'Atomic magnetic moment' in line:
            magnetic_moment.append(line.split()[-1])
            #atom_index.append(line.split()[-3])
        elif 'SPIN  1' in line:
            spin_flag = 'up'
        elif 'SPIN  2' in line:
            spin_flag = 'down'
        elif 'eigenvalues:' in line:
            eigenvalues = lines[index+1]
            if spin_flag == 'up':
                eigenvalues_up.append(eigenvalues)
            elif spin_flag == 'down':
                eigenvalues_down.append(eigenvalues)
        elif 'occupation matrix ns' in line:
            matrix_lines = len(lines[index+1].split())
            occupation_matrix = lines[index+1: index+matrix_lines+1]
            if spin_flag == 'up':
                occupation_matrix_up.append(occupation_matrix)
            elif spin_flag == 'down':
                occupation_matrix_down.append(occupation_matrix)
                
    keys = ['atom_index', 'magnetic_moment', 'eigenvalues_up', 'eigenvalues_down', 'occupation_matrix_up', 'occupation_matrix_down', 'trace_occupation']
    tuple_list = list(zip(atom_index, magnetic_moment, eigenvalues_up, eigenvalues_down, occupation_matrix_up, occupation_matrix_down, trace_occupation))
    occupations = [dict(zip(keys, values)) for values in tuple_list]

    return occupations

