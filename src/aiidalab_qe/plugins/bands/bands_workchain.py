import numpy as np
from aiida.plugins import DataFactory, WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiidalab_qe.plugins.utils import set_component_resources
from aiida import orm
from aiida.engine import ToContext, WorkChain, calcfunction

GAMMA = "\u0393"

PwBandsWorkChain = WorkflowFactory("quantumespresso.pw.bands")
ProjwfcBandsWorkChain = WorkflowFactory("wannier90_workflows.projwfcbands")


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
    for (start, end), label_start, label_end in zip(branches, labels[:-1], labels[1:]):
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

def determine_symmetry_path(structure):
    # Tolerance for checking equality
    cell_lengths = structure.cell_lengths
    cell_angles = structure.cell_angles
    tolerance = 1e-3

    # Define symmetry conditions and their corresponding types in a dictionary
    symmetry_conditions = {
        (
            math.isclose(cell_lengths[0], cell_lengths[1], abs_tol=tolerance)
            and math.isclose(cell_angles[2], 120.0, abs_tol=tolerance)
        ): "hexagonal",
        (
            math.isclose(cell_lengths[0], cell_lengths[1], abs_tol=tolerance)
            and math.isclose(cell_angles[2], 90.0, abs_tol=tolerance)
        ): "square",
        (
            math.isclose(cell_lengths[0], cell_lengths[1], abs_tol=tolerance) == False
            and math.isclose(cell_angles[2], 90.0, abs_tol=tolerance)
        ): "rectangular",
        (
            math.isclose(
                cell_lengths[1] * math.cos(math.radians(cell_angles[2])),
                cell_lengths[0] / 2,
                abs_tol=tolerance,
            )
        ): "rectangular_centered",
        (
            math.isclose(cell_lengths[0], cell_lengths[1], abs_tol=tolerance) == False
            and math.isclose(cell_angles[2], 90.0, abs_tol=tolerance) == False
        ): "oblique",
    }

    # Check for symmetry type based on conditions
    for condition, symmetry_type in symmetry_conditions.items():
        if condition:
            if symmetry_type == "rectangular_centered" or "oblique":
                cos_gamma = np.array(structure.cell[0]).dot(structure.cell[1]) / (
                    cell_lengths[0] * cell_lengths[1]
                )
                gamma = np.arccos(cos_gamma)
                eta = (1 - (cell_lengths[0] / cell_lengths[1]) * cos_gamma) / (
                    2 * np.power(np.sin(gamma), 2)
                )
                nu = 0.5 - (eta * cell_lengths[1] * cos_gamma) / cell_lengths[0]
                return generate_2d_path(symmetry_type, eta, nu)

            return generate_2d_path(symmetry_type)
        else:
            raise ValueError("Invalid symmetry type")

class BandsWorkChain(WorkChain):
    "Workchain to compute the electronic band structure"
    label = "bands"

    @classmethod
    def define(cls, spec):
        super().define(spec)
        spec.expose_inputs(PwBandsWorkChain, namespace="bands")
        spec.expose_inputs(ProjwfcBandsWorkChain, namespace="bands_projwfc")

        spec.expose_outputs(PwBandsWorkChain)
        spec.expose_outputs(ProjwfcBandsWorkChain)

        spec.outline(cls.setup, cls.run_bands, cls.results)

        spec.exit_code(400, "ERROR_WORKCHAIN_FAILED", message="The workchain bands failed.")


    @classmethod
    def get_builder_from_protocol(
        cls,
        pw_code,
        structure,
        simulation_mode,
        protocol=None,
        projwfc_code=None,
        overrides=None,
        **kwargs,
    ):
        """ Return a BandsWorkChain builder prepopulated with inputs following the specified protocol
        
        :param structure: the ``StructureData`` instance to use.
        :param pw_code: the ``Code`` instance configured for the ``quantumespresso.pw`` plugin.
        :param protocol: protocol to use, if not specified, the default will be used.
        :param projwfc_code: the ``Code`` instance configured for the ``quantumespresso.projwfc`` plugin.
        :param simulation_mode: hat type of simulation to run normal band or fat bands.
        
        """

        builder = cls.get_builder()

        if "simulation_mode" == "normal":

            args = (pw_code, structure, protocol)
            builder = PwBandsWorkChain.get_builder_from_protocol(
                *args, overrides=overrides, **kwargs
            )
        elif "simulation_mode" == "fat_bands":
            args = (pw_code, projwfc_code, structure, protocol, projwfc_code)
            builder = ProjwfcBandsWorkChain.get_builder_from_protocol(
                *args, overrides=overrides, **kwargs
            )

        if structure.pbc != (True, True, True):
            kpoints_distance = overrides["scf"]["kpoints_distance"]
            if structure.pbc == (True, False, False):
                kpoints = generate_kpath_1d(structure, kpoints_distance)
            elif structure.pbc == (True, True, False):
                kpoints = generate_kpath_2d(
                    structure, kpoints_distance, determine_symmetry_path(structure)
                )
            builder.pop("bands_kpoints_distance")
            builder.update({"bands_kpoints": kpoints})

        return builder
            

    def setup(self):
        """Define the current workchain"""
        if "bands" in self.inputs:
            self.ctx.workchain = PwBandsWorkChain
        elif "bands_projwfc" in self.inputs:
            self.ctx.workchain = ProjwfcBandsWorkChain
        else:
            self.report("No bands workchain specified")
            return self.exit_codes.ERROR_WORKCHAIN_FAILED
        
    def run_bands(self):
        """Run the bands workchain"""
        inputs = self.exposed_inputs(self.ctx.workchain)
        inputs.update(self.inputs[self.ctx.workchain.get_link_label("bands")])
        return ToContext(
            bands=self.ctx.workchain.run.get_submission_node().get_outgoing().one().node
        )

    def results(self):
        """Attach the bands results"""
        workchain = self.ctx.workchain

        if not workchain.is_finished_ok:
            self.report("Bands workchain failed")
            return self.exit_codes.ERROR_WORKCHAIN_FAILED
        else:
            self.out_many(self.ctx.bands.outputs)
            self.report("Bands workchain completed successfully")
