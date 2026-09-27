from pathlib import Path

import ipywidgets as ipw
import traitlets
from anywidget import AnyWidget


class ProgressBar(AnyWidget):
    _esm = Path(__file__).parent / "static" / "progress_bar.js"
    _css = Path(__file__).parent / "static" / "progress_bar.css"

    description = traitlets.Unicode("").tag(sync=True)
    description_layout = traitlets.Dict(default_value={}).tag(sync=True)
    value = traitlets.Float(0.0).tag(sync=True)
    bar_style = traitlets.Unicode("").tag(sync=True)
    animating = traitlets.Bool(False).tag(sync=True)

    def __init__(self, description_layout=None, **kwargs):
        if description_layout is None:
            description_layout = ipw.Layout(width="auto", flex="2 1 auto")
        elif isinstance(description_layout, dict):
            description_layout = ipw.Layout(**description_layout)

        description_layout = {
            key: value
            for key, value in description_layout.get_state().items()
            if not key.startswith("_") and value is not None
        }

        super().__init__(**kwargs)
        self.description_layout = description_layout

    @traitlets.validate("value")
    def _validate_value(self, proposal):
        if not 0 <= proposal["value"] <= 1.0:
            raise traitlets.TraitError("The value must be between 0 and 1.0.")
        return proposal["value"]
