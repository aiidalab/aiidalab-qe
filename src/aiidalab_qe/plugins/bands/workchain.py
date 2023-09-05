from aiida.plugins import WorkflowFactory

PwBandsWorkChain = WorkflowFactory("quantumespresso.pw.bands")


def get_builder(codes, structure, overrides, protocol, **kwargs):
    pw_code = codes.get("pw_code", {})
    bands_overrides = overrides.get("bands", {})
    bands = PwBandsWorkChain.get_builder_from_protocol(
        code=pw_code,
        structure=structure,
        protocol=protocol,
        overrides=bands_overrides,
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
