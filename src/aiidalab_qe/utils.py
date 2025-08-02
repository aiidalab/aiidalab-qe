import typing as t

from aiida_pseudo.groups.family import PseudoPotentialFamily
from upf_tools import UPFDict

from aiida import orm
from aiida.common.exceptions import NotExistent
from aiida.plugins import DataFactory, GroupFactory

UpfData = DataFactory("pseudo.upf")
SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")


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

RELATIVISTIC_MAPPING = {
    "SR": "scalar",
    "FR": "full",
}


def get_pseudo_by_filename(filename: str):
    return (
        orm.QueryBuilder()
        .append(
            UpfData,
            filters={"attributes.filename": filename},
        )
        .first(flat=True)
    )


def get_pseudo_by_md5(md5: str):
    return (
        orm.QueryBuilder()
        .append(
            UpfData,
            filters={"attributes.md5": md5},
        )
        .first(flat=True)
    )


def get_pseudo_family(md5: str) -> t.Optional[PseudoPotentialFamily]:
    return (
        orm.QueryBuilder()
        .append(
            UpfData,
            filters={"attributes.md5": md5},
            tag="pseudo",
        )
        .append(
            (
                SsspFamily,
                PseudoDojoFamily,
                CutoffsPseudoPotentialFamily,
            ),
            with_node="pseudo",
        )
        .first(flat=True)
    )


def get_info_from_family(pseudo_family: PseudoPotentialFamily):
    try:
        if isinstance(pseudo_family, SsspFamily):
            functional = pseudo_family.label.split("/")[2]
            relativistic = None
        else:
            functional = pseudo_family.label.split("/")[2]
            relativistic = pseudo_family.label.split("/")[3]
    except TypeError as err:
        raise ValueError(
            "Pseudo family label format is invalid. Expected format: 'family/functional/relativistic'."
        ) from err
    except Exception as err:
        raise ValueError("Invalid pseudo family label format") from err
    else:
        return functional, RELATIVISTIC_MAPPING.get(relativistic or "", relativistic)


def get_upf_dict(pseudo):
    with pseudo.as_path() as pseudo_path:
        upf_dict = UPFDict.from_upf(pseudo_path.as_posix())
    return upf_dict


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


def populate_extras(pp_uuid):
    try:
        pp_node = orm.load_node(pp_uuid)
    except NotExistent as err:
        raise ValueError(
            f"Pseudo potential with UUID {pp_uuid} does not exist"
        ) from err

    try:
        upf_dict = get_upf_dict(pp_node)
    except Exception as err:
        raise ValueError(
            f"Failed to read UPF data from pseudo potential with UUID {pp_uuid}"
        ) from err

    if pseudo_family := get_pseudo_family(pp_node.md5):
        functional, relativistic = get_info_from_family(pseudo_family)
    else:
        functional = FUNCTIONAL_MAPPING.get(
            functional := "+".join(
                upf_dict["header"].get("functional", "").lower().split()
            ),
            functional,
        )
        relativistic = upf_dict["header"].get("relativistic", None)

    pp_node.base.extras.set_many(
        {
            "functional": functional or None,
            "relativistic": relativistic,
        }
    )

    for key in ("ecutrho", "ecutwfc"):
        if key in pp_node.base.extras.all:
            pp_node.base.extras.delete(key)
