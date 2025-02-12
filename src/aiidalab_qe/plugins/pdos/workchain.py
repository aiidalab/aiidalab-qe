from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiidalab_qe.utils import (
    enable_pencil_decomposition,
    set_component_resources,
)

PdosWorkChain = WorkflowFactory("quantumespresso.pdos")
PwBandsWorkChain = WorkflowFactory("quantumespresso.pw.bands")


def check_codes(pw_code, dos_code, projwfc_code):
    """Check that the codes are installed on the same computer."""
    if (
        not any(
            [
                pw_code is None,
                dos_code is None,
                projwfc_code is None,
            ]
        )
        and len(
            {
                pw_code.computer.pk,
                dos_code.computer.pk,
                projwfc_code.computer.pk,
            }
        )
        != 1
    ):
        raise ValueError(
            "All selected codes must be installed on the same computer. This is because the "
            "PDOS calculations rely on large files that are not retrieved by AiiDA."
        )


def update_resources(builder, codes):
    set_component_resources(builder.scf.pw, codes.get("pw"))
    set_component_resources(builder.nscf.pw, codes.get("pw"))
    set_component_resources(builder.dos, codes.get("dos"))
    set_component_resources(builder.projwfc, codes.get("projwfc"))
    enable_pencil_decomposition(builder.scf.pw)
    enable_pencil_decomposition(builder.nscf.pw)
    enable_pencil_decomposition(builder.dos)

    # disable the parallelization setting for projwfc
    # npool = codes["pw"]["parallelization"]["npool"]
    # builder.projwfc.settings = orm.Dict(dict={"cmdline": ["-nk", str(npool)]})


def get_builder(codes, structure, parameters, **kwargs):
    from copy import deepcopy

    pw_code = codes.get("pw")["code"]
    dos_code = codes.get("dos")["code"]
    projwfc_code = codes.get("projwfc")["code"]
    check_codes(pw_code, dos_code, projwfc_code)
    protocol = parameters["workchain"]["protocol"]

    scf_overrides = deepcopy(parameters["advanced"])
    nscf_overrides = deepcopy(parameters["advanced"])

    # Dos Projwfc overrides
    dos_overrides = {
        "parameters": {
            "DOS": {
                "DeltaE": parameters["pdos"]["energy_grid_step"],
            }
        }
    }
    projwfc_overrides = {
        "parameters": {
            "PROJWFC": {
                "DeltaE": parameters["pdos"]["energy_grid_step"],
            }
        }
    }

    if parameters["pdos"]["use_pdos_degauss"]:
        dos_overrides["parameters"]["DOS"] = {
            "degauss": parameters["pdos"]["pdos_degauss"]
        }
        projwfc_overrides["parameters"]["PROJWFC"] = {
            "degauss": parameters["pdos"]["pdos_degauss"]
        }
        # Using Gaussian smearing with 'tetrahedra_opt' in NSCF calculations
        # causes projwfc.x to fail in producing projections.
        # This issue occurs in 3D, 2D, and 1D systems, as well as in molecules.
        # To avoid this, we use 'tetrahedra_lin' instead when 'pdos_degauss' is set.
        nscf_overrides["pw"]["parameters"]["SYSTEM"]["occupations"] = "tetrahedra_lin"

    # Update the nscf kpoints distance from the setting panel
    nscf_overrides["kpoints_distance"] = parameters["pdos"]["nscf_kpoints_distance"]

    overrides = {
        "scf": scf_overrides,
        "nscf": nscf_overrides,
        "dos": dos_overrides,
        "projwfc": projwfc_overrides,
    }
    if dos_code is not None and projwfc_code is not None:
        pdos = PdosWorkChain.get_builder_from_protocol(
            pw_code=pw_code,
            dos_code=dos_code,
            projwfc_code=projwfc_code,
            structure=structure,
            protocol=protocol,
            electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
            spin_type=SpinType(parameters["workchain"]["spin_type"]),
            initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
            overrides=overrides,
            **kwargs,
        )
        # include nbands_factor (Same as in BandsWorkChain)
        pdos["nbands_factor"] = orm.Float(
            PwBandsWorkChain.get_protocol_inputs()["nbands_factor"]
        )
        # pop the inputs that are exclueded from the expose_inputs
        pdos.pop("structure", None)
        pdos.pop("clean_workdir", None)
        # update resources
        update_resources(pdos, codes)

        if (
            scf_overrides["pw"]["parameters"]["SYSTEM"].get("tot_magnetization")
            is not None
        ):
            pdos.scf["pw"]["parameters"]["SYSTEM"].pop("starting_magnetization", None)
            pdos.nscf["pw"]["parameters"]["SYSTEM"].pop("starting_magnetization", None)
    else:
        raise ValueError("The dos_code and projwfc_code are required.")
    return pdos


def update_inputs(inputs, ctx):
    """Update the inputs using context."""
    inputs.structure = ctx.current_structure


workchain_and_builder = {
    "workchain": PdosWorkChain,
    "exclude": ("structure", "relax"),
    "get_builder": get_builder,
    "update_inputs": update_inputs,
}
