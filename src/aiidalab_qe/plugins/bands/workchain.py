from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

PwBandsWorkChain = WorkflowFactory("quantumespresso.pw.bands")


def get_builder(codes, structure, parameters, **kwargs):
    """Get a builder for the PwBandsWorkChain."""
    from copy import deepcopy

    pw_code = codes.get("pw", {})
    protocol = parameters["workchain"]["protocol"]
    #
    scf_overrides = deepcopy(parameters["advanced"])
    bands_overrides = deepcopy(parameters["advanced"])
    bands_overrides.pop("kpoints_distance", None)
    bands_overrides["pw"]["parameters"]["SYSTEM"].pop("smearing", None)
    bands_overrides["pw"]["parameters"]["SYSTEM"].pop("degauss", None)
    overrides = {
        "scf": scf_overrides,
        "bands": bands_overrides,
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
    # pop the inputs that are excluded from the expose_inputs
    bands.pop("relax")
    bands.pop("structure", None)
    bands.pop("clean_workdir", None)
    return bands


workchain_and_builder = {
    "workchain": PwBandsWorkChain,
    "exclude": ("clean_workdir", "structure", "relax"),
    "get_builder": get_builder,
}
