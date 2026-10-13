from pathlib import Path

import traitlets
from anywidget import AnyWidget


class TableWidget(AnyWidget):
    _esm = Path(__file__).parent / "static" / "table_widget.js"
    _css = Path(__file__).parent / "static" / "table_widget.css"

    data = traitlets.List().tag(sync=True)
    selected_rows = traitlets.List().tag(sync=True)
