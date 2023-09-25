from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, SpinType

PdosWorkChain = WorkflowFactory("quantumespresso.pdos")


def get_builder(codes, structure, parameters, **kwargs):
    from copy import deepcopy

    pw_code = codes.get("pw", None)
    dos_code = codes.get("dos", None)
    projwfc_code = codes.get("projwfc", None)
    protocol = parameters["workchain"]["protocol"]
    kpoints_distance = parameters["advanced"]["kpoints_distance"]

    scf_overrides = deepcopy(parameters["advanced"])
    nscf_overrides = deepcopy(parameters["advanced"])
    if protocol == "fast" and (kpoints_distance < 0.5):
        nscf_overrides["kpoints_distance"] = kpoints_distance
    if protocol == "moderate" and (kpoints_distance < 0.1):
        nscf_overrides["kpoints_distance"] = kpoints_distance
    if protocol == "precise" and (kpoints_distance < 0.05):
        nscf_overrides["kpoints_distance"] = kpoints_distance

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
