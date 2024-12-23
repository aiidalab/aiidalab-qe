from aiida import orm
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType
from aiidalab_qe.plugins.utils import set_component_resources

BandsWorkChain = WorkflowFactory("aiidalab_qe.bands_workchain")
# from .bands_workchain import BandsWorkChain


def check_codes(pw_code, projwfc_code):
    """Check that the codes are installed on the same computer."""
    if (
        not any(
            [
                pw_code is None,
                projwfc_code is None,
            ]
        )
        and len(
            {
                pw_code.computer.pk,
                projwfc_code.computer.pk,
            }
        )
        != 1
    ):
        raise ValueError(
            "All selected codes must be installed on the same computer. This is because the "
            "BandsWorkChain calculations rely on large files that are not retrieved by AiiDA."
        )


def update_resources(builder, codes):
    if "bands" in builder:
        set_component_resources(builder.bands.scf.pw, codes.get("pw"))
        set_component_resources(builder.bands.bands.pw, codes.get("pw"))
        builder.bands.scf.pw.settings = orm.Dict({"CMDLINE": ["-pd", ".true."]})
        builder.bands.bands.pw.settings = orm.Dict({"CMDLINE": ["-pd", ".true."]})
    elif "bands_projwfc" in builder:
        set_component_resources(builder.bands_projwfc.scf.pw, codes.get("pw"))
        set_component_resources(builder.bands_projwfc.bands.pw, codes.get("pw"))
        set_component_resources(
            builder.bands_projwfc.projwfc.projwfc, codes.get("projwfc_bands")
        )
        builder.bands_projwfc.scf.pw.settings = orm.Dict({"CMDLINE": ["-pd", ".true."]})
        builder.bands_projwfc.bands.pw.settings = orm.Dict(
            {"CMDLINE": ["-pd", ".true."]}
        )


def get_builder(codes, structure, parameters, **kwargs):
    """Get a builder for the BandsWorkChain."""
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

    check_codes(pw_code, codes.get("projwfc_bands")["code"])

    overrides = {
        "scf": scf_overrides,
        "bands": bands_overrides,
        "relax": relax_overrides,
    }

    if parameters["bands"]["projwfc_bands"]:
        simulation_mode = "fat_bands"
    else:
        simulation_mode = "normal"

    bands_builder = BandsWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        projwfc_code=codes.get("projwfc_bands")["code"],
        structure=structure,
        simulation_mode=simulation_mode,
        protocol=protocol,
        electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
        spin_type=SpinType(parameters["workchain"]["spin_type"]),
        initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
        overrides=overrides,
        **kwargs,
    )
    update_resources(bands_builder, codes)

    return bands_builder


def update_inputs(inputs, ctx):
    """Update the inputs using context."""
    inputs.structure = ctx.current_structure


workchain_and_builder = {
    "workchain": BandsWorkChain,
    "exclude": ("structure", "relax"),
    "get_builder": get_builder,
    "update_inputs": update_inputs,
}
