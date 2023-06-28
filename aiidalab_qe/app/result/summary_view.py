"""Workflow summary

"""
import typing

import ipywidgets as ipw

from aiidalab_qe.app.panel import ResultPanel

FUNCTIONAL_LINK_MAP = {
    "PBE": "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.77.3865",
    "PBEsol": "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.136406",
}

PSEUDO_LINK_MAP = {
    "SSSP": "https://www.materialscloud.org/discover/sssp/table/efficiency",
    "PseudoDojo": "http://www.pseudo-dojo.org/",
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


def extract_report_parameters(builder, ui_parameters) -> dict[str, typing.Any]:
    """Extract (recover) the parameters for report from the builder and ui parameters.
    There are some parameters that are not stored in the builder, but can be extracted
    directly from the widgets, such as the ``pseudo_family`` and ``relax_type``.
    """
    report = {}
    # Extract the pw calculation parameters from the builder
    # energy_cutoff is same for all pw calculations when pseudopotentials are fixed
    # as well as the smearing settings (semaring and degauss) and scf kpoints distance
    # read from the first pw calculation of relax workflow.
    # It is safe then to extract these parameters from the first pw calculation, since the
    # builder is anyway set with subworkchain inputs even it is not run which controlled by
    # the properties inputs.
    energy_cutoff_wfc = builder.relax.base["pw"]["parameters"]["SYSTEM"]["ecutwfc"]
    energy_cutoff_rho = builder.relax.base["pw"]["parameters"]["SYSTEM"]["ecutrho"]
    occupation = builder.relax.base["pw"]["parameters"]["SYSTEM"]["occupations"]
    scf_kpoints_distance = builder.relax.base.kpoints_distance.value
    report.update(
        {
            "energy_cutoff_wfc": energy_cutoff_wfc,
            "energy_cutoff_rho": energy_cutoff_rho,
            "occupation": occupation,
            "scf_kpoints_distance": scf_kpoints_distance,
        }
    )
    if occupation == "smearing":
        report["degauss"] = builder.relax.base["pw"]["parameters"]["SYSTEM"]["degauss"]
        report["smearing"] = builder.relax.base["pw"]["parameters"]["SYSTEM"][
            "smearing"
        ]
    report["bands_kpoints_distance"] = builder.bands.bands_kpoints_distance.value
    report["nscf_kpoints_distance"] = builder.pdos.nscf.kpoints_distance.value
    report["tot_charge"] = builder.relax.base["pw"]["parameters"]["SYSTEM"].get(
        "tot_charge", 0.0
    )
    # report from parameters
    # Workflow logic
    report["relax_method"] = ui_parameters["workflow"]["relax_type"]
    report["relaxed"] = ui_parameters["workflow"]["relax_type"] != "none"
    report["bands_computed"] = ui_parameters["workflow"]["properties"]["bands"]
    report["pdos_computed"] = ui_parameters["workflow"]["properties"]["pdos"]
    # Material settings
    report["material_magnetic"] = ui_parameters["basic"]["spin_type"]
    report["electronic_type"] = ui_parameters["basic"]["electronic_type"]
    # Calculation settings
    report["protocol"] = ui_parameters["basic"]["protocol"]
    # Pseudopotential settings
    pseudo_library = ui_parameters["advanced"]["pseudo_family"]
    pseudo_library_info = pseudo_library.split("/")
    if pseudo_library_info[0] == "SSSP":
        pseudo_protocol = pseudo_library_info[3]
    elif pseudo_library_info[0] == "PseudoDojo":
        pseudo_protocol = pseudo_library_info[4]
    report["pseudo_family"] = pseudo_library_info[0]
    report["pseudo_version"] = pseudo_library_info[1]
    report["pseudo_protocol"] = pseudo_protocol
    report["pseudo_library"] = pseudo_library_info[0]
    report["functional"] = pseudo_library_info[2]
    report["pseudo_link"] = PSEUDO_LINK_MAP[pseudo_library_info[0]]
    report["functional_link"] = FUNCTIONAL_LINK_MAP[report["functional"]]
    return report


def _generate_report_html(report):
    """Read from the bulider parameters and generate a html for reporting
    the inputs for the `QeAppWorkChain`.
    """
    from importlib import resources

    from jinja2 import Environment

    from aiidalab_qe.app import static

    def _fmt_yes_no(truthy):
        return "Yes" if truthy else "No"

    env = Environment()
    env.filters.update(
        {
            "fmt_yes_no": _fmt_yes_no,
        }
    )
    template = resources.read_text(static, "workflow_summary.jinja")
    style = resources.read_text(static, "style.css")

    return env.from_string(template).render(style=style, **report)


def generate_report_html(qeapp_wc):
    """Generate a html for reporting the inputs for the `QeAppWorkChain`"""
    builder = qeapp_wc.get_builder_restart()
    ui_parameters = qeapp_wc.get_extra("ui_parameters", {})
    report = dict(extract_report_parameters(builder, ui_parameters))

    return _generate_report_html(report)


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


class SummaryView(ResultPanel):
    title = "Workflow Summary"

    def _update_view(self):
        report_html = generate_report_html(self.qeapp_node)

        self.summary_view = ipw.HTML(report_html)
        self.children = [self.summary_view]
