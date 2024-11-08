"""Class to .

Authors:

    AiiDAlab Team
"""

from __future__ import annotations

import os
import typing as t

import ipywidgets as ipw
import traitlets as tl

from aiidalab_qe.common.mixins import Confirmable
from aiidalab_qe.common.mvc import Model

DEFAULT_PARAMETERS = {}


class Panel(ipw.VBox):
    """Base class for all the panels.

    The base class has a method to return the value of all the widgets in
    the panel as a dictionary. The dictionary is used to construct the
    input file for the calculation. The class also has a method to load a dictionary to set the value of the widgets in the panel.

    title: the title to be shown in the GUI
    identifier: which plugin this panel belong to.

    """

    title = "Panel"

    def __init__(self, parent=None, identifier=None, **kwargs):
        """Initialize the panel.

        :param kwargs: keyword arguments to pass to the ipw.VBox constructor.
        """
        self.parent = parent
        self.identifier = identifier or getattr(self, "identifier", "plugin")
        super().__init__(
            children=kwargs.pop("children", []),
            **kwargs,
        )


class SettingsOutline(ipw.HBox):
    title = "Outline"
    description = ""

    def __init__(self, **kwargs):
        self.include = ipw.Checkbox(
            description=self.title,
            indent=False,
            style={"description_width": "initial"},
        )

        super().__init__(
            children=[
                self.include,
                ipw.HTML(f"""
                    <div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
                        {self.description}
                    </div>
                """),
            ],
            **kwargs,
        )


class SettingsModel(Model, Confirmable):
    title = "Model"
    dependencies: list[str] = []

    include = tl.Bool(False)
    loaded_from_process = tl.Bool(False)

    _defaults = {}

    def update(self, specific=""):
        """Updates the model.

        Parameters
        ----------
        `specific` : `str`, optional
            If provided, specifies the level of update.
        """
        pass

    def get_model_state(self) -> dict:
        """Retrieves the model current state as a dictionary."""
        raise NotImplementedError

    def set_model_state(self, parameters: dict):
        """Distributes the parameters of a loaded calculation to the model."""
        raise NotImplementedError

    def reset(self):
        """Resets the model to present defaults."""
        pass


SM = t.TypeVar("SM", bound=SettingsModel)


class SettingsPanel(Panel, t.Generic[SM]):
    title = "Settings"
    description = ""

    def __init__(self, model: SM, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        self.loading_message = LoadingWidget(f"Loading {self.identifier} settings")

        super().__init__(
            children=[self.loading_message],
            **kwargs,
        )

        self._model = model

        self.rendered = False
        self.updated = False

        self.links = []

    def render(self):
        raise NotImplementedError

    def refresh(self, specific=""):
        """Refreshes the settings panel.

        Unlinks the panel's widgets. If the panel's model is included in the
        calculation, also updates the model's defaults. Resets the model to
        these defaults if there is no input structure.

        Parameters
        ----------
        `specific` : `str`, optional
            If provided, specifies the level of refresh.
        """
        self.updated = False
        self._unsubscribe()
        if self._model.include:
            self.update(specific)
        if "PYTEST_CURRENT_TEST" in os.environ:
            # Skip resetting to avoid having to inject a structure when testing
            return
        if hasattr(self._model, "input_structure") and not self._model.input_structure:
            self._reset()

    def update(self, specific=""):
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


class ResultPanel(Panel):
    """Base class for all the result panels.

    The base class has a method to load the result of the calculation.
    And a show method to display it in the panel.
    It has a update method to update the result in the panel.
    """

    title = "Result"
    # to specify which plugins (outputs) are needed for this result panel.
    workchain_labels = []

    def __init__(self, node=None, **kwargs):
        self.node = node
        super().__init__(
            children=[
                ipw.VBox(
                    [ipw.Label(f"{self.title} not available.")],
                    layout=ipw.Layout(min_height="380px"),
                )
            ],
            **kwargs,
        )

    @property
    def outputs(self):
        """Outputs of the calculation."""
        return None if self.node is None else self.node.outputs

    def _update_view(self):
        """Update the result in the panel.

        :param result: the result of the calculation.
        """
