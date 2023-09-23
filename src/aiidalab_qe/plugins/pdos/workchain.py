from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

PdosWorkChain = WorkflowFactory("quantumespresso.pdos")


def check_codes(pw_code, dos_code, projwfc_code):
    if (
        not any(
            [
                pw_code.value is None,
                dos_code.value is None,
                projwfc_code.value is None,
            ]
        )
        and len(
            set(
                (
                    pw_code.computer.pk,
                    dos_code.computer.pk,
                    projwfc_code.computer.pk,
                )
            )
        )
        != 1
    ):
        raise ValueError(
            "All selected codes must be installed on the same computer. This is because the "
            "PDOS calculations rely on large files that are not retrieved by AiiDA."
        )


def get_builder(codes, structure, parameters, **kwargs):
    from copy import deepcopy

    pw_code = codes.get("pw", None)
    dos_code = codes.get("dos", None)
    projwfc_code = codes.get("projwfc", None)
    check_codes(pw_code, dos_code, projwfc_code)
    protocol = parameters["workchain"]["protocol"]
    #
    scf_overrides = deepcopy(parameters["advanced"])
    nscf_overrides = deepcopy(parameters["advanced"])
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
            electronic_type=ElectronicType(parameters["workchain"]["electronic_type"]),
            spin_type=SpinType(parameters["workchain"]["spin_type"]),
            initial_magnetic_moments=parameters["advanced"]["initial_magnetic_moments"],
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
