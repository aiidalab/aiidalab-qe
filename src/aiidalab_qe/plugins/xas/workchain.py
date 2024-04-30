from importlib import resources

import yaml
from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

from aiidalab_qe.plugins import xas as xas_folder

XspectraCrystalWorkChain = WorkflowFactory("quantumespresso.xspectra.crystal")
PSEUDO_TOC = yaml.safe_load(resources.read_text(xas_folder, "pseudo_toc.yaml"))
pseudo_data_dict = PSEUDO_TOC["pseudos"]
xch_elements = PSEUDO_TOC["xas_xch_elements"]


def update_resources(builder, codes):
    """Update the resources for the builder."""
    builder.core.scf.pw.metadata.options.resources = {
        "num_machines": codes.get("pw")["nodes"],
        "num_mpiprocs_per_machine": codes.get("pw")["ntasks_per_node"],
        "num_cores_per_mpiproc": codes.get("pw")["cpus_per_task"],
    }
    builder.core.scf.pw.parallelization = orm.Dict(dict=codes["pw"]["parallelization"])
    builder.core.xs_prod.xspectra.metadata.options.resources = {
        "num_machines": codes.get("xspectra")["nodes"],
        "num_mpiprocs_per_machine": codes.get("xspectra")["ntasks_per_node"],
        "num_cores_per_mpiproc": codes.get("xspectra")["cpus_per_task"],
    }


def get_builder(codes, structure, parameters, **kwargs):
    from copy import deepcopy

    protocol = parameters["workchain"]["protocol"]
    xas_parameters = parameters["xas"]
    core_hole_treatments = xas_parameters["core_hole_treatments"]
    elements_list = xas_parameters["elements_list"]
    supercell_min_parameter = xas_parameters["supercell_min_parameter"]
    pseudo_labels = xas_parameters["pseudo_labels"]
    core_wfc_data_labels = xas_parameters["core_wfc_data_labels"]
    pseudos = {}
    # Convert the pseudo and core_wfc_data node labels into nodes:
    core_wfc_data = {k: orm.load_node(v) for k, v in core_wfc_data_labels.items()}
    for element in elements_list:
        pseudos[element] = {
            k: orm.load_node(v) for k, v in pseudo_labels[element].items()
        }

    # TODO should we override the cutoff_wfc, cutoff_rho by the new pseudo?
    # In principle we should, if we know what that value is, but that would
    # require testing them first...

    # (13/10/23) I'm keeping the part about molecules in for future reference,
    # but we need to establish the protocol & backend code for XAS of molecules
    # before thinking about a workflow.
    # (22/01/24) Commented out the code for molecules, just so the option doesn't
    # appear in the UI and confuse the user.
    # is_molecule_input = (
    #     True if xas_parameters.get("structure_type") == "molecule" else False
    # )

    structure_preparation_settings = {
        "supercell_min_parameter": orm.Float(supercell_min_parameter),
        # "is_molecule_input": orm.Bool(is_molecule_input),
    }
    spglib_settings = orm.Dict({"symprec": 1.0e-3})

    pw_code = codes["pw"]["code"]
    xs_code = codes["xspectra"]["code"]
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

    builder = XspectraCrystalWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        xs_code=xs_code,
        structure=structure,
        protocol=protocol,
        pseudos=pseudos,
        elements_list=elements_list,
        core_hole_treatments=core_hole_treatments,
        core_wfc_data=core_wfc_data,
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
    builder.structure_preparation_settings = structure_preparation_settings
    # update resources
    update_resources(builder, codes)

    return builder


def update_inputs(inputs, ctx):
    """Update the inputs using context."""
    inputs.structure = ctx.current_structure


workchain_and_builder = {
    "workchain": XspectraCrystalWorkChain,
    "exclude": ("structure", "relax"),
    "get_builder": get_builder,
    "update_inputs": update_inputs,
}
