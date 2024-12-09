from __future__ import annotations

from pathlib import Path

import traitlets as tl
from bs4 import BeautifulSoup

import aiidalab_qe
from aiidalab_qe.app.utils import get_entry_items


class GuideManager(tl.HasTraits):
    active_guide = tl.Unicode("none", allow_none=True)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        guides = Path(aiidalab_qe.__file__).parent.joinpath("guides").glob("*")
        self._guides = {guide.stem: guide.absolute() for guide in guides}
        self._fetch_plugin_guides()

        self.content = BeautifulSoup()

        self.observe(
            self._on_active_guide_change,
            "active_guide",
        )

    @property
    def has_guide(self):
        return self.active_guide != "none"

    def get_guides(self):
        return [*self._guides.keys()]

    def get_guide_section_by_id(self, content_id: str):
        return self.content.find(attrs={"id": content_id})

    def _on_active_guide_change(self, _):
        guide_path = self._guides.get(self.active_guide)
        html = Path(guide_path).read_text() if guide_path else ""
        self.content = BeautifulSoup(html, "html.parser")

    def _fetch_plugin_guides(self):
        entries: dict[str, Path] = get_entry_items("aiidalab_qe.properties", "guides")
        for guides in entries.values():
            for guide in guides.glob("*"):
                self._guides[guide.stem] = guide.absolute()


guide_manager = GuideManager()
