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
    """The `InAppGuide` is used to set up toggleable in-app guides.

    Attributes
    ----------
    `manager` : `GuideManager`
        A local reference to the global guide manager.
    `identifier` : `str`, optional
        If content `children` are not provided directly, the `identifier`
        is used to fetch the corresponding guide section from the guide
        currently loaded by the guide manager.

    Raises
    ------
    `ValueError`
        If neither content `children` or a guide section `identifier` are provided.
    """

    def __init__(
        self,
        children: list | None = None,
        guide_id: str = "",
        identifier: str = "",
        classes: list[str] | None = None,
        **kwargs,
    ):
        """`InAppGuide` constructor.

        Parameters
        ----------
        `children` : `list`, optional
            The content children of this guide section.
        `guide_id` : `str`, optional
            The associated guide id to be used in conjunction with content children.
            If none provided, the widget-based guide section will be shown for all
            guides.
        `identifier` : `str`, optional
            If content `children` are not provided directly, the `identifier`
            is used to fetch the corresponding guide section from the guide
            currently loaded by the guide manager.
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
            self.guide_id = guide_id
            self.children = children
            self.identifier = None
        elif identifier:
            self.guide_id = None
            self.children = []
            self.identifier = identifier
        else:
            raise ValueError("No widgets or path identifier provided")

        self.manager.observe(
            self._on_active_guide_change,
            "active_guide",
        )

        # This manual toggle call is necessary because the guide
        # may be contained in a component that was not yet rendered
        # when a guide was selected.
        self._on_active_guide_change(None)

    def _on_active_guide_change(self, _):
        self._update_contents()
        self._toggle_guide()

    def _update_contents(self):
        """Update the contents of the guide section."""
        if not self.identifier:
            return
        html = self.manager.get_guide_section_by_id(self.identifier)
        self.children = [ipw.HTML(str(html))] if html else []

    def _toggle_guide(self):
        """Toggle the visibility of the guide section."""
        self.layout.display = (
            "flex"
            if self.children
            and (
                # file-based guide section
                (self.identifier and self.manager.has_guide)
                # widget-based guide section shown for all guides
                or (not self.guide_id and self.manager.has_guide)
                # widget-based guide section shown for a specific guide
                or self.guide_id == self.manager.active_guide
            )
            else "none"
        )
