from __future__ import annotations

import ipywidgets as ipw


class InfoBox(ipw.VBox):
    """The `InfoBox` component is used to provide additional info regarding a widget or an app."""

    def __init__(self, classes: list[str] | None = None, **kwargs):
        """`InfoBox` constructor.

        Parameters
        ----------
        `classes` : `list[str]`, optional
            One or more CSS classes.
        """
        super().__init__(**kwargs)
        self.add_class("info-box")
        for custom_classes in classes or []:
            for custom_class in custom_classes.split(" "):
                if custom_class:
                    self.add_class(custom_class)


class InAppGuide(InfoBox):
    """The `InfoAppGuide` is used to set up toggleable in-app guides."""

    def __init__(
        self,
        children: list | None = None,
        identifier: str = "",
        classes: list[str] | None = None,
        **kwargs,
    ):
        """`InAppGuide` constructor.

        Parameters
        ----------
        children : `list`, optional
            The children of the guide.
        `identifier` : `str`
            The identifier used to load the guide file.
        `classes` : `list[str]`, optional
            One or more CSS classes.
        """
        from aiidalab_qe.common.guide_manager import guide_manager

        self.manager = guide_manager

        super().__init__(
            classes=[
                "in-app-guide",
                *(classes or []),
            ],
            **kwargs,
        )

        if children:
            self.children = children
        elif identifier:
            self.children = []
            self.identifier = identifier
        else:
            raise ValueError("No content or path identifier provided")

        self.manager.observe(
            self._on_active_guide_change,
            "active_guide",
        )

        # This manual toggle call is necessary because the guide
        # may be contained in a component that was not yet rendered
        # when a guide was selected.
        self._toggle_guide()

    def _on_active_guide_change(self, _):
        self._toggle_guide()

    def _toggle_guide(self):
        if hasattr(self, "identifier"):
            html = self.manager.get_guide_section_by_id(self.identifier)
            self.children = [ipw.HTML(str(html))] if html else []
        self.layout.display = "flex" if self.manager.has_guide else "none"
