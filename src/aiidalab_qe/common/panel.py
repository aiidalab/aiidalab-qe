"""Class to .

Authors:

    AiiDAlab Team
"""

from __future__ import annotations

import typing as t

import ipywidgets as ipw
import traitlets as tl

if t.TYPE_CHECKING:
    from aiidalab_qe.app.configuration.model import ConfigurationModel

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

    def reset(self):
        raise NotImplementedError


class PanelOutline(Panel):
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


class SettingPanel(Panel):
    title = "Settings"
    description = ""

    def __init__(self, config_model: ConfigurationModel, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget(f"Loading {self.identifier} settings")],
            **kwargs,
        )

        self._config_model = config_model
        self._model = config_model.get_model(self.identifier)

        self.links = []

        self.rendered = False

    def _unsubscribe(self):
        for link in self.links:
            link.unlink()
        self.links.clear()


class SettingsModel(tl.HasTraits):
    title = "Model"

    include = tl.Bool()
    confirmed = tl.Bool(False)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.observe(self._unconfirm, tl.All)

    def update(self):
        raise NotImplementedError

    def get_model_state(self):
        raise NotImplementedError

    def set_model_state(self, parameters):
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError

    def _unconfirm(self, change):
        if change["name"] != "confirmed":
            self.confirmed = False


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
        if self.node is None:
            return None

        return self.node.outputs

    def _update_view(self):
        """Update the result in the panel.

        :param result: the result of the calculation.
        """
