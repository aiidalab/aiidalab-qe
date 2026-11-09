from collections.abc import Mapping
from importlib import resources

import yaml

from aiidalab_qe.config import CONFIG_DIR_PATH

from ..setup.codes import QE_VERSION


def recursive_merge(d1, d2):
    """Merge dictionary `d2` into dictionary `d1` recursively.

    - For keys that are in both dictionaries:
        - If values are dictionaries, merge them recursively.
        - Otherwise, overwrite the value in `d1` with the value in `d2`.
    - Keys in `d2` not in `d1` are added to `d1`.

    Parameters
    ----------
    `d1` : `dict`
        The dictionary to merge into.
    `d2` : `dict`
        The dictionary to merge from.

    Returns
    -------
    `dict`
        The merged dictionary.

    Examples
    --------
    >>> d1 = {'a': 1, 'b': {'c': 2, 'd': 3}}
    >>> d2 = {'b': {'c': 4, 'e': 5}}
    >>> recursive_merge(d1, d2)
    {'a': 1, 'b': {'c': 4, 'd': 3, 'e': 5}}
    """
    for key, value in d2.items():
        if key in d1:
            if isinstance(d1[key], Mapping) and isinstance(value, Mapping):
                recursive_merge(d1[key], value)
            else:
                d1[key] = value
        else:
            d1[key] = value
    return d1


DefaultCodesDict = dict[str, dict[str, str]]


def _set_default_code_labels(default_codes: DefaultCodesDict) -> DefaultCodesDict:
    """Set the default code labels with the QE version."""
    for code_name in default_codes:
        label = default_codes[code_name]["code"]
        assert "{QE_VERSION}" in label, "Code label must contain '{QE_VERSION}'"
        default_codes[code_name]["code"] = label.format(QE_VERSION=QE_VERSION)
    return default_codes


DEFAULT_PARAMETERS = yaml.safe_load(resources.read_text(__package__, "qeapp.yaml"))
DEFAULT_PARAMETERS["codes"] = _set_default_code_labels(DEFAULT_PARAMETERS["codes"])

custom_config_file = CONFIG_DIR_PATH / "qe-config.yml"
if custom_config_file.exists():
    custom_config = yaml.safe_load(custom_config_file.read_text())
    DEFAULT_PARAMETERS = recursive_merge(DEFAULT_PARAMETERS, custom_config)
