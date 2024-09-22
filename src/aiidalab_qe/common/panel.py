"""Class to .

Authors:

    AiiDAlab Team
"""

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

    def get_panel_value(self):
        """Return the value of all the widgets in the panel as a dictionary.

        :return: a dictionary of the values of all the widgets in the panel.
        """
        return {}

    def set_panel_value(self, parameters):
        """Set the value of the widgets in the panel.

        :param parameters: a dictionary of the values of all the widgets in the panel.
        """
        for key, value in parameters.items():
            if key in self.__dict__:
                setattr(self, key, value)

    def reset(self):
        """Reset the panel to the default value."""
        self.set_panel_value(DEFAULT_PARAMETERS)

    def _update_state(self):
        """Update the state of the panel."""


class PanelModel(tl.HasTraits):
    title = "Model"

    def get_model_state(self):
        raise NotImplementedError

    def set_model_state(self, state: dict):
        raise NotImplementedError

    def reset(self):
        raise NotImplementedError


class SettingPanel(Panel):
    title = "Settings"
    description = ""

    def __init__(self, model, **kwargs):
        from aiidalab_qe.common.widgets import LoadingWidget

        super().__init__(
            children=[LoadingWidget(f"Loading {self.identifier} settings")],
            **kwargs,
        )

        self._model = model

        self.rendered = False


class OutlinePanel(Panel):
    title = "Outline"
    description = ""

    def __init__(self, **kwargs):
        # Checkbox to see if the property should be calculated
        self.run = ipw.Checkbox(
            description=self.title,
            indent=False,
            value=False,
            style={"description_width": "initial"},
        )
        self.description_html = ipw.HTML(
            f"""<div style="line-height: 140%; padding-top: 0px; padding-bottom: 5px">
            {self.description}</div>"""
        )
        self.accordion = ipw.Accordion(children=[self.description_html])
        self.accordion.selected_index = None
        super().__init__(children=[self.run, self.description_html], **kwargs)

    def get_panel_value(self):
        return {f"{self.identifier}_run": self.run.value}

    def set_panel_value(self, input_dict):
        self.run.value = input_dict.get(f"{self.identifier}_run", False)


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
