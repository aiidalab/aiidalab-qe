import os
import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.common.mvc import Model


class AdvancedSubModel(Model):
    dependencies = []

    loaded_from_process = tl.Bool(False)

    _defaults = {}

    def update(self, specific=""):
        """Updates the model.

        Parameters
        ----------
        `specific` : `str`, optional
            If provided, specifies the level of update.
        """
        raise NotImplementedError()

    def reset(self):
        """Resets the model to present defaults."""
        raise NotImplementedError()


M = t.TypeVar("M", bound=AdvancedSubModel)


class AdvancedSubSettings(ipw.VBox, t.Generic[M]):
    identifier = "sub"

    def __init__(self, model: M, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        self.loading_message = LoadingWidget(f"Loading {self.identifier} settings")

        super().__init__(
            layout={"justify_content": "space-between", **kwargs.get("layout", {})},
            children=[self.loading_message],
            **kwargs,
        )

        self._model = model
        self._model.observe(
            self._on_override_change,
            "override",
        )

        self.rendered = False
        self.updated = False

        self.links = []

    def render(self):
        raise NotImplementedError()

    def refresh(self, specific=""):
        """Refreshes the subsettings section.

        Unlinks any linked widgets and updates the model's defaults.
        Resets the model to these defaults if there is no input structure.

        Parameters
        ----------
        `specific` : `str`, optional
            If provided, specifies the level of refresh.
        """
        self.updated = False
        self._unsubscribe()
        self._update(specific)
        if "PYTEST_CURRENT_TEST" in os.environ:
            # Skip resetting to avoid having to inject a structure when testing
            return
        if hasattr(self._model, "input_structure") and not self._model.input_structure:
            self._reset()

    def _on_override_change(self, change):
        if not change["new"]:
            self._reset()

    def _update(self, specific=""):
        """Updates the model if not yet updated.

        Parameters
        ----------
        `specific` : `str`, optional
            If provided, specifies the level of update.
        """
        if self.updated:
            return
        if not self._model.loaded_from_process:
            self._model.update(specific)
        self.updated = True

    def _unsubscribe(self):
        """Unlinks any linked widgets."""
        for link in self.links:
            link.unlink()
        self.links.clear()

    def _reset(self):
        """Resets the model to present defaults."""
        self.updated = False
        self._model.reset()
