from __future__ import annotations

import ipywidgets as ipw
import traitlets as tl


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


class GuideManager(tl.HasTraits):
    active_guide = tl.Unicode("none")


guide_manager = GuideManager()


class InAppGuide(InfoBox):
    """The `InfoAppGuide` is used to set up in-app guides that may be toggle in unison."""

    def __init__(
        self,
        guide_class: str = "qe-app",
        classes: list[str] | None = None,
        **kwargs,
    ):
        """`InAppGuide` constructor.

        Parameters
        ----------
        `guide_class` : `str`, optional
            The identifier used to toggle the guide.
            The default `qe-app` identifies built-in guide sections.
        `classes` : `list[str]`, optional
            One or more CSS classes.
        """

        self.guide_class = guide_class

        super().__init__(
            classes=[
                "in-app-guide",
                *(classes or []),
            ],
            **kwargs,
        )

        guide_manager.observe(
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
        active_guide = guide_manager.active_guide
        not_generic = self.guide_class != "qe-app"
        if active_guide == "none" or (not_generic and active_guide != self.guide_class):
            self.layout.display = "none"
        else:
            self.layout.display = "flex"
