import ipywidgets as ipw
from aiida_quantumespresso.workflows.pw.bands import PwBandsWorkChain

FUNCTIONAL_LINK_MAP = {
    "PBE": "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.77.3865",
    "PBEsol": "https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.100.136406",
}

PSEUDO_LINK_MAP = {
    "SSSP": "https://www.materialscloud.org/discover/sssp/table/efficiency",
    "PseudoDojo": "http://www.pseudo-dojo.org/",
}

FUNCTIONAL_REPORT_MAP = {
    "LDA": "local density approximation (LDA)",
    "PBE": "generalized gradient approximation of Perdew-Burke-Ernzerhof (PBE)",
    "PBEsol": "the revised generalized gradient approximation of Perdew-Burke-Ernzerhof (PBE) for solids",
}

# Periodicity
PERIODICITY_MAPPING = {
    (True, True, True): "xyz",
    (True, True, False): "xy",
    (True, False, False): "x",
}


def generate_report_parameters(qeapp_wc):
    """Generate the report parameters from the ui parameters and workchain's input.

    Parameters extracted from ui parameters, directly from the widgets,
    such as the ``pseudo_family`` and ``relax_type``.

    Parameters extracted from workchain's inputs, such as the ``energy_cutoff_wfc``
    and ``energy_cutoff_rho``.

    Return a dictionary of the parameters.
    """
    from aiida.orm.utils.serialize import deserialize_unsafe

    ui_parameters = qeapp_wc.base.extras.get("ui_parameters", {})
    if isinstance(ui_parameters, str):
        ui_parameters = deserialize_unsafe(ui_parameters)
    # Construct the report parameters needed for the report
    # drop support for old ui parameters
    if "workchain" not in ui_parameters:
        return {}
    report = {
        "relaxed": ui_parameters["workchain"]["relax_type"],
        "relax_method": ui_parameters["workchain"]["relax_type"],
        "electronic_type": ui_parameters["workchain"]["electronic_type"],
        "material_magnetic": ui_parameters["workchain"]["spin_type"],
        "protocol": ui_parameters["workchain"]["protocol"],
        "initial_magnetic_moments": ui_parameters["advanced"][
            "initial_magnetic_moments"
        ],
        "properties": ui_parameters["workchain"]["properties"],
    }
    #
    report.update(
        {
            "bands_computed": "bands" in ui_parameters["workchain"]["properties"],
            "pdos_computed": "pdos" in ui_parameters["workchain"]["properties"],
        }
    )
    # update pseudo family information to report
    pseudo_family = ui_parameters["advanced"].get("pseudo_family")
    pseudo_family_info = pseudo_family.split("/")
    pseudo_library = pseudo_family_info[0]
    functional = pseudo_family_info[2]
    if pseudo_library == "SSSP":
        pseudo_protocol = pseudo_family_info[3]
    elif pseudo_library == "PseudoDojo":
        pseudo_protocol = pseudo_family_info[4]
    report.update(
        {
            "pseudo_family": pseudo_family,
            "pseudo_library": pseudo_library,
            "pseudo_version": pseudo_family_info[1],
            "functional": functional,
            "pseudo_protocol": pseudo_protocol,
            "pseudo_link": PSEUDO_LINK_MAP[pseudo_library],
            "functional_link": FUNCTIONAL_LINK_MAP[functional],
        }
    )
    # Extract the pw calculation parameters from the workchain's inputs
    # energy_cutoff is same for all pw calculations when pseudopotentials are fixed
    # as well as the smearing settings (semaring and degauss) and scf kpoints distance
    # read from the first pw calculation of relax workflow.
    # It is safe then to extract these parameters from the first pw calculation, since the
    # builder is anyway set with subworkchain inputs even it is not run which controlled by
    # the properties inputs.
    pw_parameters = qeapp_wc.inputs.relax.base.pw.parameters.get_dict()
    energy_cutoff_wfc = pw_parameters["SYSTEM"]["ecutwfc"]
    energy_cutoff_rho = pw_parameters["SYSTEM"]["ecutrho"]
    occupation = pw_parameters["SYSTEM"]["occupations"]
    scf_kpoints_distance = qeapp_wc.inputs.relax.base.kpoints_distance.value
    report.update(
        {
            "energy_cutoff_wfc": energy_cutoff_wfc,
            "energy_cutoff_rho": energy_cutoff_rho,
            "occupation_type": occupation,
            "scf_kpoints_distance": scf_kpoints_distance,
        }
    )
    if occupation == "smearing":
        report["degauss"] = pw_parameters["SYSTEM"]["degauss"]
        report["smearing"] = pw_parameters["SYSTEM"]["smearing"]
    report["tot_charge"] = pw_parameters["SYSTEM"].get("tot_charge", 0.0)
    report["periodicity"] = PERIODICITY_MAPPING.get(
        qeapp_wc.inputs.structure.pbc, "xyz"
    )

    # DFT+U
    hubbard_dict = ui_parameters["advanced"].pop("hubbard_parameters", None)
    if hubbard_dict:
        hubbard_parameters = hubbard_dict["hubbard_u"]
        report["hubbard_u"] = hubbard_parameters
    # hard code bands and pdos
    if "bands" in qeapp_wc.inputs:
        report["bands_kpoints_distance"] = PwBandsWorkChain.get_protocol_inputs(
            report["protocol"]
        )["bands_kpoints_distance"]

    if "pdos" in qeapp_wc.inputs:
        report[
            "nscf_kpoints_distance"
        ] = qeapp_wc.inputs.pdos.nscf.kpoints_distance.value
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
    report = {key: value for key, value in report.items() if value is not None}

    return env.from_string(template).render(style=style, **report)


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


class SummaryView(ipw.VBox):
    def __init__(self, wc_node, **kwargs):
        self.report = generate_report_parameters(wc_node)
        self.report_html = _generate_report_html(self.report)

        self.summary_view = ipw.HTML(self.report_html)
        super().__init__(
            children=[self.summary_view],
            **kwargs,
        )
