from collections.abc import Mapping
from importlib import resources
from pathlib import Path

import yaml

from aiidalab_qe.app import parameters


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


DEFAULT_PARAMETERS = yaml.safe_load(resources.read_text(parameters, "qeapp.yaml"))

CONFIG_DIR_PATH = Path.home() / ".aiidalab" / "quantum-espresso"

custom_config_file = CONFIG_DIR_PATH / "qe-config.yml"
if custom_config_file.exists():
    custom_config = yaml.safe_load(custom_config_file.read_text())
    DEFAULT_PARAMETERS = recursive_merge(DEFAULT_PARAMETERS, custom_config)
