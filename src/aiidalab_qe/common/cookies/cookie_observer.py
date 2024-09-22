from __future__ import annotations

import ipywidgets as ipw
from traitlets import Bool, Unicode


@ipw.register
class Cookie(ipw.DOMWidget):
    _view_name = Unicode("CookieObserverView").tag(sync=True)
    _view_module = Unicode("cookie_observer").tag(sync=True)
    _view_module_version = Unicode("0.1.0").tag(sync=True)
    name = Unicode("").tag(sync=True)
    value = Bool(None, allow_none=True).tag(sync=True)

    def __init__(self, name: str, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.name = name
        self.layout.display = "none"

    @property
    def exists(self) -> bool:
        return bool(self.value)
