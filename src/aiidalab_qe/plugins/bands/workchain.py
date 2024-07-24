import numpy as np
from aiida.plugins import DataFactory, WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

from aiidalab_qe.plugins.utils import set_component_resources

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
    return max(
        2, int(np.round(distance / bands_kpoints_distance))
    )  # at least two points for each segment, including both endpoints explicitly


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
    path = selected_paths[kpath_2d]["path"]
    labels = selected_paths[kpath_2d]["labels"]
    branches = zip(path[:-1], path[1:])

    all_kpoints = []  # List to hold all k-points
    label_map = []  # List to hold labels and their corresponding k-point indices

    # Calculate the number of points per branch and generate the kpoints
    index_offset = 0  # Start index for each segment
    for (start, end), label_start, _label_end in zip(branches, labels[:-1], labels[1:]):
        num_points_per_branch = points_per_branch(
            start, end, reciprocal_cell, bands_kpoints_distance
        )
        # Exclude endpoint except for the last segment to prevent duplication
        points = np.linspace(start, end, num=num_points_per_branch, endpoint=False)
        all_kpoints.extend(points)
        label_map.append(
            (index_offset, label_start)
        )  # Label for the start of the segment
        index_offset += len(points)

    # Include the last point and its label
    all_kpoints.append(path[-1])
    label_map.append((index_offset, labels[-1]))  # Label for the last point

    # Set the kpoints and their labels in KpointsData
    kpoints.set_kpoints(all_kpoints)
    kpoints.labels = label_map

    return kpoints


def update_resources(builder, codes):
    set_component_resources(builder.scf.pw, codes.get("pw"))
    set_component_resources(builder.bands.pw, codes.get("pw"))


def get_builder(codes, structure, parameters, **kwargs):
    """Get a builder for the PwBandsWorkChain."""
    from copy import deepcopy

    pw_code = codes.get("pw")["code"]
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
    # update resources
    update_resources(bands, codes)

    if scf_overrides["pw"]["parameters"]["SYSTEM"].get("tot_magnetization") is not None:
        bands.scf["pw"]["parameters"]["SYSTEM"].pop("starting_magnetization", None)
        bands.bands["pw"]["parameters"]["SYSTEM"].pop("starting_magnetization", None)

    return bands


def update_inputs(inputs, ctx):
    """Update the inputs using context."""
    inputs.structure = ctx.current_structure
    inputs.scf.pw.parameters = inputs.scf.pw.parameters.get_dict()
    if ctx.current_number_of_bands:
        inputs.scf.pw.parameters.setdefault("SYSTEM", {}).setdefault(
            "nbnd", ctx.current_number_of_bands
        )


workchain_and_builder = {
    "workchain": PwBandsWorkChain,
    "exclude": ("structure", "relax"),
    "get_builder": get_builder,
    "update_inputs": update_inputs,
}
