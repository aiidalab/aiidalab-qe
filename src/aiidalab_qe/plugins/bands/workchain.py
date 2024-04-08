import numpy as np
from aiida.plugins import DataFactory, WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

GAMMA = "\u0393"

PwBandsWorkChain = WorkflowFactory("quantumespresso.pw.bands")
KpointsData = DataFactory("core.array.kpoints")


def points_per_branch(vector_a, vector_b, reciprocal_cell, bands_kpoints_distance):
    """function to calculate the number of points per branch depending on the kpoints_distance and the reciprocal cell"""
    scaled_vector_a = np.array(vector_a)
    scaled_vector_b = np.array(vector_b)
    reciprocal_vector_a = scaled_vector_a.dot(reciprocal_cell)
    reciprocal_vector_b = scaled_vector_b.dot(reciprocal_cell)
    distance = np.linalg.norm(reciprocal_vector_a - reciprocal_vector_b)
    return round(distance / bands_kpoints_distance)


def calculate_bands_kpoints_distance(kpoints_distance):
    """function to calculate the bands_kpoints_distance depending on the kpoints_distance"""
    if kpoints_distance >= 0.5:
        return 0.1
    elif 0.15 < kpoints_distance < 0.5:
        return 0.025
    else:
        return 0.015


def generate_kpath_1d(structure, kpoints_distance):
    """Return a kpoints object for one dimensional systems (from Gamma to X)
    The number of kpoints is calculated based on the kpoints_distance (as in the PwBandsWorkChain protocol)
    """
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    reciprocal_cell = kpoints.reciprocal_cell
    bands_kpoints_distance = calculate_bands_kpoints_distance(kpoints_distance)

    # Number of points per branch
    num_points_per_branch = points_per_branch(
        [0.0, 0.0, 0.0],
        [0.5, 0.0, 0.0],
        reciprocal_cell,
        bands_kpoints_distance,
    )
    # Generate the kpoints
    points = np.linspace(
        start=[0.0, 0.0, 0.0],
        stop=[0.5, 0.0, 0.0],
        endpoint=True,
        num=num_points_per_branch,
    )
    kpoints.set_kpoints(points.tolist())
    kpoints.labels = [[0, GAMMA], [len(points) - 1, "X"]]
    return kpoints


def generate_kpath_2d(structure, kpoints_distance, kpath_2d):
    """Return a kpoints object for two dimensional systems based on the selected 2D symmetry path
    The number of kpoints is calculated based on the kpoints_distance (as in the PwBandsWorkChain protocol)
    The 2D symmetry paths are defined as in The Journal of Physical Chemistry Letters 2022 13 (50), 11581-11594 (https://pubs.acs.org/doi/10.1021/acs.jpclett.2c02972)
    """
    kpoints = KpointsData()
    kpoints.set_cell_from_structure(structure)
    reciprocal_cell = kpoints.reciprocal_cell
    bands_kpoints_distance = calculate_bands_kpoints_distance(kpoints_distance)

    # dictionary with the 2D symmetry paths
    selected_paths = {
        "hexagonal": {
            "path": [
                [0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0],
                [0.33333, 0.33333, 0.0],
                [1.0, 0.0, 0.0],
            ],
            "labels": [GAMMA, "M", "K", GAMMA],
        },
        "square": {
            "path": [
                [0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [1.0, 0.0, 0.0],
            ],
            "labels": [GAMMA, "X", "M", GAMMA],
        },
        "rectangular": {
            "path": [
                [0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0],
                [0.5, 0.5, 0.0],
                [0.0, 0.5, 0.0],
                [1.0, 0.0, 0.0],
            ],
            "labels": [GAMMA, "X", "S", "Y", GAMMA],
        },
    }
    # if the selected path is centered_rectangular or oblique, the path is calculated based on the reciprocal cell
    if kpath_2d in ["centered_rectangular", "oblique"]:
        a1 = reciprocal_cell[0]
        a2 = reciprocal_cell[1]
        norm_a1 = np.linalg.norm(a1)
        norm_a2 = np.linalg.norm(a2)
        cos_gamma = (
            a1.dot(a2) / (norm_a1 * norm_a2)
        )  # Angle between a1 and a2 # like in https://pubs.acs.org/doi/10.1021/acs.jpclett.2c02972
        gamma = np.arccos(cos_gamma)
        eta = (1 - (norm_a1 / norm_a2) * cos_gamma) / (2 * np.power(np.sin(gamma), 2))
        nu = 0.5 - (eta * norm_a2 * cos_gamma) / norm_a1
        selected_paths["centered_rectangular"] = {
            "path": [
                [0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0],
                [1 - eta, nu, 0],
                [0.5, 0.5, 0.0],
                [eta, 1 - nu, 0.0],
                [1.0, 0.0, 0.0],
            ],
            "labels": [GAMMA, "X", "H_1", "C", "H", GAMMA],
        }
        selected_paths["oblique"] = {
            "path": [
                [0.0, 0.0, 0.0],
                [0.5, 0.0, 0.0],
                [1 - eta, nu, 0],
                [0.5, 0.5, 0.0],
                [eta, 1 - nu, 0.0],
                [0.0, 0.5, 0.0],
                [1.0, 0.0, 0.0],
            ],
            "labels": [GAMMA, "X", "H_1", "C", "H", "Y", GAMMA],
        }
    points_branch = []
    num_per_branch = []
    path = selected_paths[kpath_2d]["path"]
    labels = selected_paths[kpath_2d]["labels"]
    branches = zip(
        path[:-1], path[1:]
    )  # zip the path with the next point in the path to define the branches

    # Calculate the number of points per branch and generate the kpoints
    for branch in branches:
        num_points_per_branch = points_per_branch(
            branch[0], branch[1], reciprocal_cell, bands_kpoints_distance
        )
        if branch[1] == [1.0, 0.0, 0.0]:
            points = np.linspace(
                start=branch[0],
                stop=branch[1],
                endpoint=True,
                num=num_points_per_branch,
            )
        else:
            points = np.linspace(
                start=branch[0], stop=branch[1], num=num_points_per_branch
            )
        points_branch.append(points.tolist())
        num_per_branch.append(num_points_per_branch)

    # Generate the kpoints as single list and add the labels
    list_kpoints = [item for sublist in points_branch for item in sublist]
    kpoints.set_kpoints(list_kpoints)
    kpoints.labels = [
        [index, labels[index]]
        if index == 0
        else [list_kpoints.index(value, 1), labels[index]]
        for index, value in enumerate(path)
    ]
    return kpoints


def get_builder(codes, structure, parameters, **kwargs):
    """Get a builder for the PwBandsWorkChain."""
    from copy import deepcopy

    pw_code = codes.get("pw")
    protocol = parameters["workchain"]["protocol"]
    scf_overrides = deepcopy(parameters["advanced"])
    relax_overrides = {
        "base": deepcopy(parameters["advanced"]),
        "base_final_scf": deepcopy(parameters["advanced"]),
    }
    bands_overrides = deepcopy(parameters["advanced"])
    bands_overrides.pop("kpoints_distance", None)
    bands_overrides["pw"]["parameters"]["SYSTEM"].pop("smearing", None)
    bands_overrides["pw"]["parameters"]["SYSTEM"].pop("degauss", None)
    overrides = {
        "scf": scf_overrides,
        "bands": bands_overrides,
        "relax": relax_overrides,
    }
    bands = PwBandsWorkChain.get_builder_from_protocol(
        code=pw_code,
        structure=structure,
        protocol=protocol,
        electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
        spin_type=SpinType(parameters["workchain"]["spin_type"]),
        initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
        overrides=overrides,
        **kwargs,
    )

    if structure.pbc != (True, True, True):
        kpoints_distance = parameters["advanced"]["kpoints_distance"]
        if structure.pbc == (True, False, False):
            kpoints = generate_kpath_1d(structure, kpoints_distance)
        elif structure.pbc == (True, True, False):
            kpoints = generate_kpath_2d(
                structure, kpoints_distance, parameters["bands"]["kpath_2d"]
            )
        bands.pop("bands_kpoints_distance")
        bands.update({"bands_kpoints": kpoints})

    # pop the inputs that are excluded from the expose_inputs
    bands.pop("relax")
    bands.pop("structure", None)
    bands.pop("clean_workdir", None)

    if scf_overrides["pw"]["parameters"]["SYSTEM"].get("tot_magnetization") is not None:
        bands.scf["pw"]["parameters"]["SYSTEM"].pop("starting_magnetization", None)
        bands.bands["pw"]["parameters"]["SYSTEM"].pop("starting_magnetization", None)

    return bands


workchain_and_builder = {
    "workchain": PwBandsWorkChain,
    "exclude": ("structure", "relax"),
    "get_builder": get_builder,
}
