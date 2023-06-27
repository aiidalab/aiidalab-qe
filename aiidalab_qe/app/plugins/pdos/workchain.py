from aiida.orm import load_code
from aiida.plugins import WorkflowFactory

PdosWorkChain = WorkflowFactory("quantumespresso.pdos")


def get_builder(codes, structure, parameters):
    protocol = parameters["basic"].pop("protocol", "fast")
    pw_code = load_code(codes.get("pw_code"))
    dos_code = load_code(codes.get("dos_code"))
    projwfc_code = load_code(codes.get("projwfc_code"))
    overrides = {
        "scf": parameters["advanced"],
        "nscf": parameters["advanced"],
    }
    parameters = parameters["basic"]
    builder = PdosWorkChain.get_builder_from_protocol(
        pw_code=pw_code,
        dos_code=dos_code,
        projwfc_code=projwfc_code,
        structure=structure,
        protocol=protocol,
        overrides=overrides,
        **parameters,
    )
    # pop the inputs that are exclueded from the expose_inputs
    builder.pop("structure", None)
    builder.pop("clean_workdir", None)
    return builder


workchain_and_builder = [PdosWorkChain, get_builder]
