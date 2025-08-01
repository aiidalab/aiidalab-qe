import contextlib

from aiida_pseudo.groups.family import PseudoPotentialFamily
from upf_tools import UPFDict

from aiida import orm


def set_component_resources(component, code_info):
    """Set the resources for a given component based on the code info."""
    if code_info:  # Ensure code_info is not None or empty (# XXX: ? from jyu, need to pop a warning to plugin developer or what?)
        code: orm.Code = code_info["code"]
        if code.computer.scheduler_type == "hyperqueue":
            component.metadata.options.resources = {
                "num_cpus": code_info["nodes"]
                * code_info["ntasks_per_node"]
                * code_info["cpus_per_task"]
            }
        else:
            # XXX: jyu should properly deal with None type of scheduler_type which can be "core.direct" (will be replaced by hyperqueue) and "core.slurm" ...
            component.metadata.options.resources = {
                "num_machines": code_info["nodes"],
                "num_mpiprocs_per_machine": code_info["ntasks_per_node"],
                "num_cores_per_mpiproc": code_info["cpus_per_task"],
            }

        component.metadata.options["max_wallclock_seconds"] = code_info[
            "max_wallclock_seconds"
        ]
        if "parallelization" in code_info:
            component.parallelization = orm.Dict(dict=code_info["parallelization"])


def enable_pencil_decomposition(component):
    """Enable the pencil decomposition for the given component."""

    component.settings = orm.Dict({"CMDLINE": ["-pd", ".true."]})


def fetch_pseudo_family_by_label(label) -> PseudoPotentialFamily:
    """Fetch the pseudo family by label."""
    return orm.Group.collection.get(label=label)  # type: ignore


def shallow_copy_nested_dict(d):
    """Recursively copies only the dictionary structure but keeps value references."""
    if isinstance(d, dict):
        return {key: shallow_copy_nested_dict(value) for key, value in d.items()}
    return d


FUNCTIONAL_MAPPING = {
    "pbe": "PBE",
    "sla+pw+pbx+pbc": "PBE",
    "pbesol": "PBEsol",
    "sla+pw+psx+psc": "PBEsol",
}


def get_upf_dict(pseudo):
    with pseudo.as_path() as pseudo_path:
        upf_dict = UPFDict.from_upf(pseudo_path.as_posix())
    return upf_dict


def get_pseudo_info(pp_uuid):
    functional = "unknown"
    relativistic = "unknown"
    ecutwfc = 0.0
    ecutrho = 0.0
    with contextlib.suppress(Exception):
        pp_node = orm.load_node(pp_uuid)
        keys = ("functional", "relativistic", "ecutwfc", "ecutrho")
        if all(key in pp_node.base.extras for key in keys):
            functional = pp_node.base.extras.get("functional")
            relativistic = pp_node.base.extras.get("relativistic")
            ecutwfc = pp_node.base.extras.get("ecutwfc")
            ecutrho = pp_node.base.extras.get("ecutrho")
        else:
            upf_dict = get_upf_dict(pp_node)
            functional = FUNCTIONAL_MAPPING.get(
                functional := "+".join(
                    upf_dict["header"]["functional"].lower().split()
                ),
                functional,
            )
            relativistic = upf_dict["header"]["relativistic"]
            ecutwfc = float(int(upf_dict["header"]["wfc_cutoff"]))
            ecutrho = float(int(upf_dict["header"]["rho_cutoff"]))
    return {
        "functional": functional,
        "relativistic": relativistic,
        "ecutwfc": ecutwfc,
        "ecutrho": ecutrho,
    }
