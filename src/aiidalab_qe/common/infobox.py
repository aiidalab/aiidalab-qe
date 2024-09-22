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
    """The `InfoAppGuide` is used to set up in-app guides that may be toggle in unison."""

    def __init__(self, guide_id: str, classes: list[str] | None = None, **kwargs):
        """`InAppGuide` constructor.

        Parameters
        ----------
        `guide_id` : `str`
            The unique identifier for the guide.
        `classes` : `list[str]`, optional
            One or more CSS classes.
        """
        classes = ["in-app-guide", *(classes or []), guide_id]
        super().__init__(classes=classes, **kwargs)


class FirstVisitBox(InfoBox):
    """The `FirstVisitBox` is used to display a message to first time users."""

    def __init__(
        self,
        message: str = "Welcome to the app",
        classes: list[str] | None = None,
        **kwargs,
    ):
        """`FirstVisitBox` constructor.

        Parameters
        ----------
        `message` : `str`, optional
            The message to display.
        `classes` : `list[str]`, optional
            One or more CSS classes.
        """
        from aiidalab_qe.common.cookies import CookieManager

        self.close_button = ipw.Button(
            icon="times",
            tooltip="Close",
        )

        self.message_box = ipw.HBox(
            children=[
                ipw.HTML(message),
                self.close_button,
            ],
        )
        self.message_box.add_class("message-box")

        self.undo_button = ipw.Button(
            icon="undo",
            tooltip="Undo",
            description="undo",
        )

        self.closing_message = ipw.HBox(
            children=[
                ipw.HTML("This message will be hidden on your next visit"),
                self.undo_button,
            ],
        )
        self.closing_message.add_class("closing-message")

        self.cookie_manager = CookieManager()
        self.cookie_name = "aiidalab-qe-user"
        self.cookie_manager.add_cookie(
            cookie_name=self.cookie_name,
            callback=self._toggle_display,
        )

        super().__init__(
            classes=classes,
            children=[self.message_box],
            **kwargs,
        )

        self.add_class("first-visit-box")

        self._set_event_listeners()

    def _toggle_display(self, _=None):
        """Toggle the first time user info box."""
        if self.cookie_manager.has_cookie(self.cookie_name):
            self.layout.display = "none"
        else:
            self.layout.display = "flex"

    def _on_close(self, _=None):
        """Display closing message and add cookie to browser."""
        self.children = [self.closing_message]
        self.cookie_manager.set_cookie(self.cookie_name)

    def _on_undo(self, _=None):
        """Display first time user message and remove cookie from browser."""
        self.children = [self.message_box]
        self.cookie_manager.remove_cookie(self.cookie_name)

    def _set_event_listeners(self):
        """Set event listeners."""
        self.close_button.on_click(self._on_close)
        self.undo_button.on_click(self._on_undo)
