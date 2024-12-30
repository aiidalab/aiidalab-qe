from aiida_pseudo.groups.family import PseudoPotentialFamily

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
