from pathlib import Path

from IPython.display import Javascript, display

display(Javascript((Path(__file__).parent / "cookie_observer.js").read_text()))

from .cookie_manager import CookieManager

__all__ = [
    "CookieManager",
]
