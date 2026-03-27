from pathlib import Path

import traitlets as tl
from anywidget import AnyWidget


class LinkButton(AnyWidget):
    _esm = Path(__file__).parent / "link_button.js"
    _css = Path(__file__).parent / "link_button.css"

    icon = tl.Unicode("").tag(sync=True)
    description = tl.Unicode("Open").tag(sync=True)
    link = tl.Unicode("").tag(sync=True)
    tooltip = tl.Unicode("").tag(sync=True)
    target = tl.Unicode("_blank").tag(sync=True)
    class_ = tl.Unicode("").tag(sync=True)
    style_ = tl.Unicode("").tag(sync=True)
    disabled = tl.Bool(False).tag(sync=True)
    prevent_default = tl.Bool(False).tag(sync=True)
    clicks = tl.Int(0).tag(sync=True)

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def on_click(self, callback, remove=False):
        if remove:
            self.unobserve(callback, "clicks")
        else:
            self.observe(callback, "clicks")
