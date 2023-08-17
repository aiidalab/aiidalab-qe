# trigger registration of the viewer widget:
from .node_view import CalcJobNodeViewerWidget  # noqa: F401
from .process import WorkChainSelector

__all__ = [
    "WorkChainSelector",
]
