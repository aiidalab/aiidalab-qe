from __future__ import annotations

from pathlib import Path

import traitlets as tl
from bs4 import BeautifulSoup, PageElement

import aiidalab_qe
from aiidalab_qe.app.utils import get_entry_items


class GuideManager(tl.HasTraits):
    """A global guide manager that loads and manages guide sections."""

    active_guide = tl.Unicode("none", allow_none=True)

    def __init__(self, *args, **kwargs):
        """`GuideManager` constructor."""

        super().__init__(*args, **kwargs)
        guides = Path(aiidalab_qe.__file__).parent.joinpath("guides").glob("*")
        self._guides = {
            "general": {
                guide.stem.split("_", maxsplit=1)[1]: guide.absolute()
                for guide in sorted(guides, key=lambda x: x.stem.split("_")[0])
            }
        }

        self._fetch_plugin_guides()

        self.content = BeautifulSoup()

        self.observe(
            self._on_active_guide_change,
            "active_guide",
        )

    @property
    def has_guide(self) -> bool:
        return self.active_guide != "none"

    def get_guide_categories(self) -> list[str]:
        """Return a list of available guides.

        Returns
        -------
        `list[str]`
            A list of the names of available guides.
        """
        return [*self._guides.keys()]

    def get_guides(self, identifier: str) -> list[str]:
        """Return a list of available sub-guides.

        Returns
        -------
        `list[str]`
            A list of the names of available sub-guides.
        """
        return [*self._guides[identifier].keys()] if identifier != "none" else []

    def get_guide_section_by_id(self, content_id: str) -> PageElement | None:
        """Return a guide section by its HTML `id` attribute.

        Parameters
        ----------
        `content_id` : `str`
            The HTML `id` attribute of the guide section.

        Returns
        -------
        `PageElement` | `None`
            The guide section or `None` if not found.
        """
        return self.content.find(attrs={"id": content_id})  # type: ignore

    def _on_active_guide_change(self, _):
        """Load the contents of the active guide."""
        if self.active_guide == "none":
            self.content = BeautifulSoup()
            return
        category, guide = self.active_guide.split("/")
        guide_path = self._guides[category][guide]
        html = Path(guide_path).read_text() if guide_path else ""
        self.content = BeautifulSoup(html, "html.parser")

    def _fetch_plugin_guides(self):
        """Fetch guides from plugins."""
        entries: dict = get_entry_items("aiidalab_qe.properties", "guides")
        for identifier, guides in entries.items():
            path = Path(guides)
            if identifier not in self._guides:
                self._guides[identifier] = {}
            for guide in sorted(path.glob("*"), key=lambda x: x.stem.split("_")[0]):
                stem = guide.stem.split("_", maxsplit=1)[1]
                self._guides[identifier][stem] = guide.absolute()


guide_manager = GuideManager()
