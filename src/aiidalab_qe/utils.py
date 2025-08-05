import typing as t

from aiida import orm
from aiida.common.exceptions import NotExistent


def generate_alert(alert_type: str, message: str, class_: str = "", style_: str = ""):
    return f"""
        <div class="alert alert-{alert_type} warning-box {class_}" style="{style_}">
            {message}
        </div>
    """


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


def shallow_copy_nested_dict(d):
    """Recursively copies only the dictionary structure but keeps value references."""
    if isinstance(d, dict):
        return {key: shallow_copy_nested_dict(value) for key, value in d.items()}
    return d


def get_pseudo_info(pp_uuid) -> dict[str, t.Any]:
    try:
        pp_node = orm.load_node(pp_uuid)
    except NotExistent as err:
        raise ValueError(
            f"Pseudo potential with UUID {pp_uuid} does not exist"
        ) from err

    return {
        "functional": pp_node.base.extras.get("functional", None),
        "relativistic": pp_node.base.extras.get("relativistic", None),
    }
