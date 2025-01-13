from __future__ import annotations

import json
import re

import numpy as np
from pymatgen.core.periodic_table import Element

from aiida.common.extendeddicts import AttributeDict
from aiida.orm import ProjectionData, WorkChainNode

# Constants for HTML tags
HTML_TAGS = {
    "s": "s",
    "pz": "p<sub>z</sub>",
    "px": "p<sub>x</sub>",
    "py": "p<sub>y</sub>",
    "dz2": "d<sub>z<sup>2</sup></sub>",
    "dxy": "d<sub>xy</sub>",
    "dxz": "d<sub>xz</sub>",
    "dyz": "d<sub>yz</sub>",
    "dx2-y2": "d<sub>x<sup>2</sup>-y<sup>2</sup></sub>",
    "fz3": "f<sub>z<sup>3</sup></sub>",
    "fxz2": "f<sub>xz<sup>2</sup></sub>",
    "fyz2": "f<sub>yz<sup>2</sup></sub>",
    "fxyz": "f<sub>xzy</sub>",
    "fx(x2-3y2)": "f<sub>x(x<sup>2</sup>-3y<sup>2</sup>)</sub>",
    "fy(3x2-y2)": "f<sub>y(3x<sup>2</sup>-y<sup>2</sup>)</sub>",
    "fy(x2-z2)": "f<sub>y(x<sup>2</sup>-z<sup>2</sup>)</sub>",
    0.5: "<sup>+1</sup>/<sub>2</sub>",
    -0.5: "<sup>-1</sup>/<sub>2</sub>",
    1.5: "<sup>+3</sup>/<sub>2</sub>",
    -1.5: "<sup>-3</sup>/<sub>2</sub>",
    2.5: "<sup>+5</sup>/<sub>2</sub>",
    -2.5: "<sup>-5</sup>/<sub>2</sub>",
    " @ ": "<br>",
    " @ [": "-[",
    "l": "<i>l</i>",
    "m_j": "m<sub>j</sub>",
}


def extract_pdos_output(node: WorkChainNode) -> AttributeDict | None:
    """Extract the PDOS output node from the given node.

    Parameters
    ----------
    `node`: `WorkChainNode`
        The node to extract the PDOS output from.

    Returns
    -------
    `AttributeDict | None`
        The PDOS output node, if available.
    """
    if not node:
        return
    if node.process_label == "QeAppWorkChain" and "pdos" in node.outputs:
        return node.outputs.pdos
    if "dos" in node.outputs and "projwfc" in node.outputs:
        items = {key: getattr(node.outputs, key) for key in node.outputs}
        return AttributeDict(items)


def extract_bands_output(node: WorkChainNode) -> AttributeDict | None:
    """Extract the bands output node from the given node.

    Parameters
    ----------
    `node`: `WorkChainNode`
        The node to extract the bands output from.

    Returns
    -------
    `AttributeDict | None`
        The bands output node, if available.
    """
    if not node:
        return
    if node.process_label == "QeAppWorkChain" and "bands" in node.outputs:
        outputs = node.outputs.bands
    else:
        outputs = node.outputs
    return (
        outputs.bands
        if "bands" in outputs
        else outputs.bands_projwfc
        if "bands_projwfc" in outputs
        else None
    )


def get_bands_data(outputs, fermi_energy=None):
    if "band_structure" not in outputs:
        return None

    bands_data = outputs.band_structure._get_bandplot_data(
        cartesian=True, prettify_format=None, join_symbol=None, get_segments=True
    )
    # The fermi energy from band calculation is not robust.
    if "fermi_energy_up" in outputs.band_parameters:
        bands_data["fermi_energy_up"] = outputs.band_parameters["fermi_energy_up"]
        bands_data["fermi_energy_down"] = outputs.band_parameters["fermi_energy_down"]
    else:
        bands_data["fermi_energy"] = (
            fermi_energy
            if fermi_energy is not None
            else outputs.band_parameters["fermi_energy"]
        )

    bands_data["pathlabels"] = _get_bands_labeling(bands_data)

    return bands_data


def get_bands_projections_data(
    outputs,
    bands_data,
    group_tag,
    plot_tag,
    selected_atoms,
    bands_width,
):
    if "projwfc" not in outputs:
        return None

    projections = []
    bands_projection = {}

    if "projections" in outputs.projwfc:
        projections.append(
            _projections_curated_options(
                outputs.projwfc.projections,
                spin_type="none",
                group_tag=group_tag,
                plot_tag=plot_tag,
                selected_atoms=selected_atoms,
                projections_pdos="projections",
            )
        )
    else:
        for spin_proj, spin_type in zip(
            [
                outputs.projwfc.projections_up,
                outputs.projwfc.projections_down,
            ],
            ["up", "down"],
        ):
            projections.append(
                _projections_curated_options(
                    spin_proj,
                    spin_type=spin_type,
                    group_tag=group_tag,
                    plot_tag=plot_tag,
                    selected_atoms=selected_atoms,
                    projections_pdos="projections",
                )
            )

    bands_projection["projected_bands"] = _prepare_projections_to_plot(
        bands_data, projections, bands_width
    )
    if plot_tag != "total":
        band_parameters: dict = outputs.band_parameters.get_dict()
        if not band_parameters.get("spin_orbit_calculation"):
            bands_projection["projected_bands"] = _update_pdos_labels(
                bands_projection["projected_bands"]
            )
    return bands_projection["projected_bands"]


def get_pdos_data(pdos, group_tag, plot_tag, selected_atoms):
    dos = []

    if "output_dos" not in pdos.dos:
        return None

    _, energy_dos, _ = pdos.dos.output_dos.get_x()
    tdos_values = {f"{n}": v for n, v, _ in pdos.dos.output_dos.get_y()}

    if "projections" in pdos.projwfc:
        # Total DOS
        tdos = {
            "label": "Total DOS",
            "x": energy_dos.tolist(),
            "y": tdos_values.get("dos").tolist(),
            "borderColor": "#8A8A8A",  # dark gray
            "backgroundColor": "#999999",  # light gray
            "backgroundAlpha": "40%",
            "lineStyle": "solid",
        }
        dos.append(tdos)
        dos += _projections_curated_options(
            pdos.projwfc.projections,
            spin_type="none",
            group_tag=group_tag,
            plot_tag=plot_tag,
            selected_atoms=selected_atoms,
        )
    else:
        # Total DOS (↑) and Total DOS (↓)
        tdos_up = {
            "label": "Total DOS (↑)",
            "x": energy_dos.tolist(),
            "y": tdos_values.get("dos_spin_up").tolist(),
            "borderColor": "#8A8A8A",  # dark gray
            "backgroundColor": "#999999",  # light gray
            "backgroundAlpha": "40%",
            "lineStyle": "solid",
        }
        tdos_down = {
            "label": "Total DOS (↓)",
            "x": energy_dos.tolist(),
            "y": (-tdos_values.get("dos_spin_down")).tolist(),
            "borderColor": "#8A8A8A",  # dark gray
            "backgroundColor": "#999999",  # light gray
            "backgroundAlpha": "40%",
            "lineStyle": "dash",
        }
        dos += [tdos_up, tdos_down]

        # Spin-up (↑) and Spin-down (↓)
        dos += _projections_curated_options(
            pdos.projwfc.projections_up,
            spin_type="up",
            group_tag=group_tag,
            plot_tag=plot_tag,
            selected_atoms=selected_atoms,
        )
        dos += _projections_curated_options(
            pdos.projwfc.projections_down,
            spin_type="down",
            line_style="dash",
            group_tag=group_tag,
            plot_tag=plot_tag,
            selected_atoms=selected_atoms,
        )

    data_dict = {
        "dos": dos,
    }
    if "fermi_energy_up" in pdos.nscf.output_parameters:
        data_dict["fermi_energy_up"] = pdos.nscf.output_parameters["fermi_energy_up"]
        data_dict["fermi_energy_down"] = pdos.nscf.output_parameters[
            "fermi_energy_down"
        ]
    else:
        data_dict["fermi_energy"] = pdos.nscf.output_parameters["fermi_energy"]

    # Updata labels if plot_tag is different than total and SOC is false
    if plot_tag != "total":
        output_parameters: dict = pdos.nscf.output_parameters.get_dict()
        if not output_parameters.get("spin_orbit_calculation", False):
            data_dict = _update_pdos_labels(data_dict)

    return json.loads(json.dumps(data_dict))


def prepare_combined_plotly_traces(x_to_conc, y_to_conc):
    """Combine multiple lines into a single trace.

    The rows of y are concatenated with a np.nan column as a separator. Moreover,
    the x values are ajduced to match the shape of the concatenated y values. These
    transfomred arrays, representing multiple datasets/lines, can be plotted in a single trace.
    """
    if y_to_conc.ndim != 2:
        raise ValueError("y must be a 2D array")

    y_dim0 = y_to_conc.shape[0]

    # Add a np.nan column as a separator
    y_transf = np.hstack(
        [
            y_to_conc,
            np.full((y_dim0, 1), np.nan),
        ]
    ).flatten()

    # Same logic for the x axis
    x_transf = x_to_conc.reshape(1, -1) * np.ones(y_dim0).reshape(-1, 1)
    x_transf = np.hstack([x_transf, np.full((y_dim0, 1), np.nan)]).flatten()

    return x_transf, y_transf


def find_max_up_and_down(x_data_up, y_data_up, x_data_down, y_data_down, x_min, x_max):
    """
    Function to find the maximum positive value and the most negative value.

    Parameters:
    - x_data_up: List of x values for the positive part.
    - y_data_up: List of y values for the positive part.
    - x_data_down: List of x values for the negative part.
    - y_data_down: List of y values for the negative part.
    - x_min: Minimum x value for the range.
    - x_max: Maximum x value for the range.

    Returns:
    - most_negative_down: Most negative value found in the down part.
    - max_up: Maximum value found in the up part.
    """
    max_up = _find_extreme_in_range(x_data_up, y_data_up, x_min, x_max, is_max=True)
    most_negative_down = _find_extreme_in_range(
        x_data_down, y_data_down, x_min, x_max, is_max=False, initial_value=0
    )

    return most_negative_down, max_up


def find_max_in_range(x_data, y_data, x_min, x_max):
    """
    Function to find the maximum value in a specified range.

    Parameters:
    - x_data: List of x values.
    - y_data: List of y values.
    - x_min: Minimum x value for the range.
    - x_max: Maximum x value for the range.

    Returns:
    - Maximum value found in the range, or None if no valid values are found.
    """
    return _find_extreme_in_range(x_data, y_data, x_min, x_max, is_max=True)


def _get_grouping_key(
    group_tag,
    plot_tag,
    atom_position,
    kind_name,
    orbital_name_plotly,
    orbital_angular_momentum,
):
    """Generates the grouping key based on group_tag and plot_tag."""

    key_formats = {
        ("atoms", "total"): r"{var1}-{var}",
        ("kinds", "total"): r"{var1}",
        ("atoms", "orbital"): r"{var1}-{var2}<br>{var}",
        ("kinds", "orbital"): r"{var1}-{var2}",
        ("atoms", "angular_momentum"): r"{var1}-{var3}<br>{var}",
        ("kinds", "angular_momentum"): r"{var1}-{var3}",
    }

    key = key_formats.get((group_tag, plot_tag))
    if key is not None:
        return key.format(
            var=atom_position,
            var1=kind_name,
            var2=orbital_name_plotly,
            var3=orbital_angular_momentum,
        )
    else:
        return None


def _curate_orbitals(orbital):
    """Curate and transform the orbital data into the desired format."""
    orbital_data = orbital.get_orbital_dict()
    kind_name = orbital_data["kind_name"]
    atom_position = [round(i, 2) for i in orbital_data["position"]]
    radial_node = orbital_data["radial_nodes"]

    try:
        orbital_name = orbital.get_name_from_quantum_numbers(
            orbital_data["angular_momentum"], orbital_data["magnetic_number"]
        ).lower()
        orbital_name_plotly = (
            f"r{radial_node} {HTML_TAGS.get(orbital_name, orbital_name)}"
        )
        orbital_angular_momentum = f"r{radial_node} {orbital_name[0]}"

    except AttributeError:
        # Set quanutum numbers
        qn_j = orbital_data["total_angular_momentum"]
        qn_l = orbital_data["angular_momentum"]
        qn_m_j = orbital_data["magnetic_number"]
        orbital_name = f"j {qn_j} l {qn_l} m_j{qn_m_j}"
        orbital_name_plotly = f"j={HTML_TAGS.get(qn_j, qn_j)} <i>l</i>={qn_l} m<sub>j</sub>={HTML_TAGS.get(qn_m_j, qn_m_j)}"
        orbital_angular_momentum = f"<i>l</i>={qn_l} "

    return orbital_name_plotly, orbital_angular_momentum, kind_name, atom_position


def _prepare_projections_to_plot(bands_data, projections, bands_width):
    """Prepare the projected bands to be plotted.

    This function transforms the projected bands into a format that can be plotted
    in a single trace. To use the fill option `toself`,
    a band needs to be concatenated with its mirror image, first.
    """
    projected_bands = []
    for spin in [0, 1]:
        # In case of non-spin-polarized calculations, the spin index is only 0
        if spin not in bands_data["band_type_idx"]:
            continue

        x_bands = bands_data["x"]
        # New shape: (number of bands, number of kpoints)
        y_bands = bands_data["y"][:, bands_data["band_type_idx"] == spin].T

        for proj in projections[spin]:
            # Create the upper and lower boundary of the fat bands based on the orbital projections
            y_bands_proj_upper = y_bands + bands_width / 2 * proj["projections"].T
            y_bands_proj_lower = y_bands - bands_width / 2 * proj["projections"].T
            # As mentioned above, the bands need to be concatenated with their mirror image
            # to create the filled areas properly
            y_bands_mirror = np.hstack(
                [y_bands_proj_upper, y_bands_proj_lower[:, ::-1]]
            )
            # Same logic for the energy axis
            x_bands_mirror = np.concatenate([x_bands, x_bands[::-1]]).reshape(1, -1)
            x_bands_comb, y_bands_proj_comb = prepare_combined_plotly_traces(
                x_bands_mirror, y_bands_mirror
            )

            projected_bands.append(
                {
                    "x": x_bands_comb.tolist(),
                    "y": y_bands_proj_comb.tolist(),
                    "label": proj["label"],
                    "color": proj["color"],
                }
            )
    return projected_bands


def _projections_curated_options(
    projections: ProjectionData,
    group_tag,
    plot_tag,
    selected_atoms,
    projections_pdos="pdos",
    spin_type="none",
    line_style="solid",
):
    """Extract and curate the projections.

    This function can be used to extract the PDOS or the projections data.
    """
    _proj_pdos = {}
    list_positions = []

    # Constants for spin types
    SPIN_LABELS = {"up": "(↑)", "down": "(↓)", "none": ""}
    SIGN_MULT_FACTOR = {"up": 1, "down": -1, "none": 1}

    if projections_pdos == "pdos":
        proj_data = projections.get_pdos()
    elif projections_pdos == "projections":
        proj_data = projections.get_projections()
    else:
        raise ValueError(f"Invalid value for `projections_pdos`: {projections_pdos}")

    for orb_proj in proj_data:
        if projections_pdos == "pdos":
            orbital, proj_pdos, energy = orb_proj
        elif projections_pdos == "projections":
            orbital, proj_pdos = orb_proj
            energy = None

        (
            orbital_name_plotly,
            orbital_angular_momentum,
            kind_name,
            atom_position,
        ) = _curate_orbitals(orbital)

        if atom_position not in list_positions:
            list_positions.append(atom_position)

        key = _get_grouping_key(
            group_tag,
            plot_tag,
            atom_position,
            kind_name,
            orbital_name_plotly,
            orbital_angular_momentum,
        )
        if not selected_atoms:
            if key:
                _proj_pdos.setdefault(key, [energy, 0])[1] += proj_pdos

        else:
            try:
                index = list_positions.index(atom_position)
                if index in selected_atoms:
                    if key:
                        _proj_pdos.setdefault(key, [energy, 0])[1] += proj_pdos

            except ValueError:
                pass

    curated_proj = []
    for label, (energy, proj_pdos) in _proj_pdos.items():
        label += SPIN_LABELS[spin_type]  # noqa: PLW2901
        if projections_pdos == "pdos":
            orbital_proj_pdos = {
                "label": label,
                "x": energy.tolist(),
                "y": (SIGN_MULT_FACTOR[spin_type] * proj_pdos).tolist(),
                "borderColor": _cmap(label),
                "lineStyle": line_style,
            }
        else:
            orbital_proj_pdos = {
                "label": label,
                "projections": proj_pdos,
                "color": _cmap(label),
            }
        curated_proj.append(orbital_proj_pdos)

    return curated_proj


def _get_bands_labeling(bandsdata: dict) -> list:
    """Function to return two lists containing the labels and values (kpoint) for plotting.
    params:
    - bandsdata: dictionary from `get_bands_projections_data` function
    output:  update bandsdata with a new key "pathlabels" including (list of str), label_values (list of float)
    """
    UNICODE_SYMBOL = {
        "GAMMA": "\u0393",
        "DELTA": "\u0394",
        "LAMBDA": "\u039b",
        "SIGMA": "\u03a3",
        "EPSILON": "\u0395",
    }
    paths = bandsdata.get("paths")
    labels = []
    for path in paths:  # Remove duplicates
        label_a = [path["from"], path["x"][0]]
        label_b = [path["to"], path["x"][-1]]
        if label_a not in labels:
            labels.append(label_a)
        if label_b not in labels:
            labels.append(label_b)

    clean_labels = []  # Format
    for i in labels:
        if clean_labels:
            if (i not in clean_labels) and (clean_labels[-1][-1] == i[1]):
                clean_labels[-1][0] = clean_labels[-1][0] + "|" + i[0]
            else:
                clean_labels.append(i)
        else:
            clean_labels.append(i)

    path_labels = [label[0] for label in clean_labels]
    for i, label in enumerate(path_labels):
        path_labels[i] = re.sub(
            r"([A-Z]+)", lambda x: UNICODE_SYMBOL.get(x.group(), x.group()), label
        )
    path_values = [label[1] for label in clean_labels]
    return [path_labels, path_values]


def _cmap(label: str) -> str:
    """Return RGB string of color for given pseudo info
    Hardcoded at the momment.
    """
    import random

    # if a unknow type generate random color based on ascii sum
    ascn = sum([ord(c) for c in label])
    random.seed(ascn)

    return f"#{random.randint(0, 0xFFFFFF):06x}"


def _find_extreme_in_range(
    x_data, y_data, x_min, x_max, is_max=True, initial_value=float("-inf")
):
    """
    General function to find the extreme value (max or min) in a given range.

    Parameters:
    - x_data: List of x values.
    - y_data: List of y values.
    - x_min: Minimum x value for the range.
    - x_max: Maximum x value for the range.
    - is_max: Boolean to determine whether to find the maximum or minimum value.
    - initial_value: Initial value for extreme (default is -inf for max and 0 for min).

    Returns:
    - Extreme value found in the range, or None if no valid values are found.
    """
    extreme_value = initial_value

    for x, y in zip(x_data, y_data):
        if x_min <= x <= x_max:
            if (is_max and y > extreme_value) or (not is_max and y < extreme_value):
                extreme_value = y

    return extreme_value if extreme_value != initial_value else None


def _get_labels_radial_nodes(pdos_dict):
    """
    Extracts the original labels from the PDOS data and constructs an orbital dictionary.

    Args:
        pdos_dict (dict): Dictionary containing PDOS data with 'dos' key representing orbital information.

    Returns:
        tuple:
            - original_labels (list): List of strings representing the original orbital labels.
            - orbital_dict (dict): A nested dictionary mapping atom kinds to orbital types and their corresponding radial nodes.
    """
    original_labels = []
    orbital_dict = {}

    label_data_list = pdos_dict["dos"] if "dos" in pdos_dict else pdos_dict
    for label_data in label_data_list:
        # for label_data in pdos_dict["dos"]:
        label_str = label_data["label"]
        original_labels.append(label_str)

        parts = label_str.split("-")
        if len(parts) < 2:
            continue  # Skip invalid or non-orbital labels

        atom = parts[0]  # Atom type (e.g., 'Fe1')
        radial_orbital = parts[1].split()  # Splits 'r# orbital' (e.g., 'r0 s')

        if len(radial_orbital) < 2:
            continue  # Malformed label

        radial_node = int(radial_orbital[0][1:])  # Extract radial index from 'r#'
        orbital = radial_orbital[1][0]  # Orbital type ('s', 'p', 'd', 'f')

        # Populate orbital_dict with atoms, orbitals, and radial nodes
        orbital_dict.setdefault(atom, {}).setdefault(orbital, set()).add(radial_node)

    return original_labels, orbital_dict


def _assign_orbital_labels(orbital_dict):
    """
    Assigns orbital labels to atoms based on their radial nodes and electronic structure.

    Args:
        orbital_dict (dict): A nested dictionary mapping atom kinds to orbital types and their corresponding radial nodes.

    Returns:
        dict: A dictionary mapping atoms and orbitals to their assigned radial labels.
    """
    result = {}

    for atom_with_number, orbitals in orbital_dict.items():
        # Extract element name (remove numeric suffixes)
        atom = re.sub(r"\d+", "", atom_with_number)
        element = Element(atom)
        electronic_structure = list(reversed(element.full_electronic_structure))

        orbital_assignment = {orb: {} for orb in ["s", "p", "d", "f"]}

        # Map orbitals from electronic structure
        orbital_map = {
            "s": [
                f"{n}{orbital}"
                for n, orbital, _ in electronic_structure
                if orbital == "s"
            ],
            "p": [
                f"{n}{orbital}"
                for n, orbital, _ in electronic_structure
                if orbital == "p"
            ],
            "d": [
                f"{n}{orbital}"
                for n, orbital, _ in electronic_structure
                if orbital == "d"
            ],
            "f": [
                f"{n}{orbital}"
                for n, orbital, _ in electronic_structure
                if orbital == "f"
            ],
        }

        # Assign radial nodes to orbitals in reverse order
        for orb_type in ["s", "p", "d", "f"]:
            if orb_type in orbitals:
                sorted_indices = sorted(orbitals[orb_type], reverse=True)
                for idx, radial_node in enumerate(sorted_indices):
                    if radial_node < len(orbital_map[orb_type]):
                        orbital_assignment[orb_type][idx] = orbital_map[orb_type][
                            radial_node
                        ][0]

        # Clean up empty orbital assignments
        result[atom_with_number] = {
            orb: val for orb, val in orbital_assignment.items() if val
        }

    return result


def _get_new_pdos_labels(input_list, orbital_dict):
    output_list = []

    for item in input_list:
        # Check if the label contains a '-' to proceed with splitting
        if "-" in item:
            before_dash, after_dash = item.split("-", 1)

            # Split the part after the dash into words to isolate the radial node (r#)
            parts = after_dash.split()

            if parts[0].startswith("r"):
                radial_index = int(parts[0][1:])  # Extract the number after 'r'

                # Check if the first element after removing the radial part corresponds to an orbital
                orbital = parts[1]

                # If the atom and orbital type exist in the orbital_dict, map the radial node
                if (
                    before_dash in orbital_dict
                    and orbital[0] in orbital_dict[before_dash]
                ):
                    if radial_index in orbital_dict[before_dash][orbital[0]]:
                        # Get the mapped radial value
                        new_radial_value = orbital_dict[before_dash][orbital[0]][
                            radial_index
                        ]

                        # Rebuild the string, removing the space before the orbital
                        after_dash = after_dash.replace(
                            f"r{radial_index}", new_radial_value, 1
                        )
                        after_dash = after_dash.replace(
                            " ", "", 1
                        )  # Remove the space after the radial value
                        new_item = f"{before_dash}-{after_dash}"
                    else:
                        new_item = (
                            item  # If radial index not found, use the original item
                        )
                else:
                    new_item = item  # If no match in orbital_dict, use original label
            else:
                new_item = item  # In case there's no valid 'r#' part
        else:
            new_item = item  # If no dash, use the original item

        output_list.append(new_item)

    return output_list


def _update_pdos_labels(pdos_data):
    """
    Updates PDOS labels by assigning correct radial nodes to orbitals based on their electronic structure.

    Args:
        pdos_data (dict): PDOS data structure containing 'dos' key with orbital information.

    Returns:
        tuple:
            - pdos_data (dict): Updated PDOS data with correct orbital labels.
    """
    original_labels, orbital_dict = _get_labels_radial_nodes(pdos_data)
    orbital_assignment = _assign_orbital_labels(orbital_dict)
    updated_labels = _get_new_pdos_labels(original_labels, orbital_assignment)

    label_data_list = pdos_data["dos"] if "dos" in pdos_data else pdos_data
    # Update labels directly using zip
    for label_data, label in zip(label_data_list, updated_labels):
        label_data["label"] = label

    return pdos_data


def hex_to_rgba(hex_code, alpha=1):
    # Remove the '#' if it's included
    hex_code = hex_code.lstrip("#")

    # Convert hex to RGB values
    r = int(hex_code[0:2], 16)
    g = int(hex_code[2:4], 16)
    b = int(hex_code[4:6], 16)

    # Return the rgba color
    return f"rgba({r}, {g}, {b}, {alpha})"


def rgba_to_hex(color_str):
    """
    Converts an RGBA color string to a HEX color string.

    Parameters:
        color_str (str): A color string in the format 'rgba(r, g, b, a)'.

    Returns:
        str: The HEX color string in the format '#RRGGBB' (ignoring alpha) or '#RRGGBBAA' (if needed).
    """
    if color_str.startswith("rgba"):
        # Extract RGBA values from the string
        rgba_values = color_str[5:-1].split(",")
        rgba_values = [
            float(value.strip()) for value in rgba_values
        ]  # Convert to floats

        # Map RGBA values to integers (alpha scaled to 255)
        r, g, b = map(int, rgba_values[:3])
        # a = int(rgba_values[3] * 255)  # Scale alpha to 0-255

        # Convert to HEX including alpha (optional)
        return f"#{r:02x}{g:02x}{b:02x}"

    elif color_str.startswith("rgb"):
        # Handle RGB format without alpha
        rgb_values = color_str[4:-1].split(",")
        rgb_values = [int(value.strip()) for value in rgb_values]  # Convert to integers

        # Map RGB values and ignore alpha
        r, g, b = rgb_values[:3]
        return f"#{r:02x}{g:02x}{b:02x}"

    # Return unchanged if already in hex or invalid format
    return color_str


def replace_html_tags(input_string, html_tags):
    """
    Replaces HTML parts in the input string with their corresponding keys from the HTML_TAGS dictionary.

    Args:
        input_string (str): The string potentially containing HTML tags.
        html_tags (dict): Dictionary mapping keys to HTML tag replacements.

    Returns:
        str: The input string with HTML tags replaced by their keys.
    """
    for key, value in html_tags.items():
        input_string = input_string.replace(value, str(key))  # Ensure key is a string
    return input_string
