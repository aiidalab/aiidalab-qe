from aiida import orm


def set_component_resources(component, code_info):
    """Set the resources for a given component based on the code info."""
    if code_info:  # Ensure code_info is not None or empty
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
