from aiida.plugins import WorkflowFactory

PdosWorkChain = WorkflowFactory("quantumespresso.pdos")


def get_builder(codes, structure, overrides, protocol, **kwargs):
    pw_code = codes.get("pw_code", None)
    dos_code = codes.get("dos_code", None)
    projwfc_code = codes.get("projwfc_code", None)
    from copy import deepcopy

    pw_code = codes.get("pw_code", {})
    scf_overrides = deepcopy(overrides)
    nscf_overrides = deepcopy(overrides)
    nscf_overrides.pop("kpoints_distance", None)
    nscf_overrides["pw"]["parameters"]["SYSTEM"].pop("smearing", None)
    nscf_overrides["pw"]["parameters"]["SYSTEM"].pop("degauss", None)
    overrides = {
        "scf": scf_overrides,
        "nscf": nscf_overrides,
    }
    if dos_code is not None and projwfc_code is not None:
        pdos = PdosWorkChain.get_builder_from_protocol(
            pw_code=pw_code,
            dos_code=dos_code,
            projwfc_code=projwfc_code,
            structure=structure,
            protocol=protocol,
            overrides=overrides,
            **kwargs,
        )
        # pop the inputs that are exclueded from the expose_inputs
        pdos.pop("structure", None)
        pdos.pop("clean_workdir", None)
    else:
        raise ValueError("The dos_code and projwfc_code are required.")
    return pdos


workchain_and_builder = {
    "workchain": PdosWorkChain,
    "exclude": ("clean_workdir", "structure", "relax"),
    "get_builder": get_builder,
}
