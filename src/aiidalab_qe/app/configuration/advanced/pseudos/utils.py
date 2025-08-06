import typing as t

from aiida_pseudo.groups.family import PseudoPotentialFamily
from upf_tools import UPFDict

from aiida import orm
from aiida.plugins import DataFactory, GroupFactory

UpfData = DataFactory("pseudo.upf")
SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")

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


def get_pseudo_family_by_content(md5: str) -> t.Optional[PseudoPotentialFamily]:
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


def get_pseudo_family_by_label(label) -> PseudoPotentialFamily:
    return t.cast(PseudoPotentialFamily, orm.Group.collection.get(label=label))


def get_upf_dict(pseudo):
    with pseudo.as_path() as pseudo_path:
        upf_dict = UPFDict.from_upf(pseudo_path.as_posix())
    return upf_dict
