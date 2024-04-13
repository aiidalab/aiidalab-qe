from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

PdosWorkChain = WorkflowFactory("quantumespresso.pdos")


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

    pw_code = codes.get("pw")
    dos_code = codes.get("dos")
    projwfc_code = codes.get("projwfc")
    check_codes(pw_code, dos_code, projwfc_code)
    protocol = parameters["workchain"]["protocol"]

    scf_overrides = deepcopy(parameters["advanced"])
    nscf_overrides = deepcopy(parameters["advanced"])

    # Update the nscf kpoints distance from the setting panel
    nscf_overrides["kpoints_distance"] = parameters["pdos"]["nscf_kpoints_distance"]

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

        if (
            scf_overrides["pw"]["parameters"]["SYSTEM"].get("tot_magnetization")
            is not None
        ):
            pdos.scf["pw"]["parameters"]["SYSTEM"].pop("starting_magnetization", None)
            pdos.nscf["pw"]["parameters"]["SYSTEM"].pop("starting_magnetization", None)
    else:
        raise ValueError("The dos_code and projwfc_code are required.")
    return pdos


workchain_and_builder = {
    "workchain": PdosWorkChain,
    "exclude": ("structure", "relax"),
    "get_builder": get_builder,
}
