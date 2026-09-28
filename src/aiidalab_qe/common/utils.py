import typing as t

import ipywidgets as ipw


def normalize_layout(layout: ipw.Layout | dict | None) -> dict[str, t.Any]:
    """Normalize an ipywidgets Layout or dictionary into a dictionary.

    :param layout: An optional ipywidgets Layout object or dictionary.
    :return: A dictionary representing the layout.
    """
    if layout is None:
        return {}

    if isinstance(layout, dict):
        return layout.copy()

    return {
        key: value
        for key, value in layout.get_state().items()
        if not key.startswith("_") and value is not None
    }
