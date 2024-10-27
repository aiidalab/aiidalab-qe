"""Class to .

Authors:

    AiiDAlab Team
"""

from __future__ import annotations

import os

import ipywidgets as ipw
import traitlets as tl

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


class SettingsModel(tl.HasTraits):
    title = "Model"
    dependencies: list[str] = []

    include = tl.Bool()
    confirmed = tl.Bool(False)

    _defaults = {}

    def __init__(self, include=False, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.observe(self._unconfirm, tl.All)
        self.include = include

    @property
    def has_pbc(self):
        if hasattr(self, "input_structure"):
            return self.input_structure is None or any(self.input_structure.pbc)
        return False

    def update(self):
        pass

    def get_model_state(self) -> dict:
        raise NotImplementedError

    def set_model_state(self, parameters: dict):
        raise NotImplementedError

    def reset(self):
        pass

    def _unconfirm(self, change):
        if change["name"] != "confirmed":
            self.confirmed = False

    def _update_defaults(self):
        raise NotImplementedError


class SettingsPanel(Panel):
    title = "Settings"
    description = ""

    def __init__(self, model: SettingsModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        self.loading_message = LoadingWidget(f"Loading {self.identifier} settings")

        super().__init__(
            children=[self.loading_message],
            **kwargs,
        )

        self._model = model

        self.links = []

        self.rendered = False
        self.updated = False

    def render(self):
        raise NotImplementedError

    def _unsubscribe(self):
        for link in self.links:
            link.unlink()
        self.links.clear()

        if not self._model.include:
            return
        self.updated = False
        self._unsubscribe()
        self._update()
        if "PYTEST_CURRENT_TEST" in os.environ:
            # Skip resetting to avoid having to inject a structure when testing
            return
        if not self._model.input_structure:
            self._model.reset()

    def _update(self):
        if self.updated:
            return
        self._model.update()
        self.updated = True


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
