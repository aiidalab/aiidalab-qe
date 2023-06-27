from aiida.orm import load_code
from aiida.plugins import WorkflowFactory
from aiida_quantumespresso.common.types import ElectronicType, RelaxType, SpinType

PwRelaxWorkChain = WorkflowFactory("quantumespresso.pw.relax")


def get_builder(codes, structure, parameters):
    protocol = parameters["basic"].pop("protocol", "fast")
    pw_code = load_code(codes.get("pw_code"))
    # TODO check the overrides
    relax_overrides = {
        "scf": parameters["advanced"],
        "bands": parameters["advanced"],
    }
    #
    relax_type = RelaxType(parameters["workflow"]["relax_type"])
    parameters["basic"]["electronic_type"] = ElectronicType(
        parameters["basic"]["electronic_type"]
    )
    parameters["basic"]["spin_type"] = SpinType(parameters["basic"]["spin_type"])
    relax_parameters = parameters["basic"]
    builder = PwRelaxWorkChain.get_builder_from_protocol(
        code=pw_code,
        structure=structure,
        protocol=protocol,
        overrides=relax_overrides,
        relax_type=relax_type,
        **relax_parameters,
    )
    # pop the inputs that are excluded from the expose_inputs
    builder.pop("structure", None)
    builder.pop("clean_workdir", None)
    builder.pop("base_final_scf", None)  # never run a final scf
    return builder


workchain_and_builder = [PwRelaxWorkChain, get_builder]
