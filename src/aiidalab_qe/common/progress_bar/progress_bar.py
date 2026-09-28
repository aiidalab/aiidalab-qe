import typing as t
from pathlib import Path

import ipywidgets as ipw
import traitlets
from anywidget import AnyWidget

from aiidalab_qe.common.utils import normalize_layout


class ProgressBar(AnyWidget):
    _esm = Path(__file__).parent / "static" / "progress_bar.js"
    _css = Path(__file__).parent / "static" / "progress_bar.css"

    description = traitlets.Unicode("").tag(sync=True)
    description_layout = traitlets.Dict(default_value={}).tag(sync=True)
    value = traitlets.Float(0.0).tag(sync=True)
    bar_style = traitlets.Unicode("").tag(sync=True)
    animating = traitlets.Bool(False).tag(sync=True)

    def __init__(
        self,
        description_layout: ipw.Layout | dict | None = None,
        **kwargs: t.Any,
    ) -> None:
        description_layout = normalize_layout(description_layout)
        description_layout.setdefault("width", "auto")
        description_layout.setdefault("flex", "1 1 auto")

        super().__init__(**kwargs)
        self.description_layout = description_layout

    @traitlets.validate("value")
    def _validate_value(self, proposal):
        if not 0 <= proposal["value"] <= 1.0:
            raise traitlets.TraitError("The value must be between 0 and 1.0.")
        return proposal["value"]
