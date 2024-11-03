from __future__ import annotations

import io

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.plugins import DataFactory, GroupFactory
from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_widgets_base.utils import StatusHTML

from ..subsettings import AdvancedSubSettings
from .model import PseudosModel

UpfData = DataFactory("pseudo.upf")
SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")


class PseudoSettings(AdvancedSubSettings):
    identifier = "pseudos"

    def __init__(self, model: PseudosModel, **kwargs):
        super().__init__(model, **kwargs)

        self._model.observe(
            self._on_input_structure_change,
            "input_structure",
        )
        self._model.observe(
            self._on_protocol_change,
            "protocol",
        )
        self._model.observe(
            self._on_spin_orbit_change,
            "spin_orbit",
        )
        self._model.observe(
            self._on_family_parameters_change,
            ["library", "functional"],
        )
        self._model.observe(
            self._on_family_change,
            "family",
        )

        ipw.dlink(
            (self._model, "cutoffs"),
            (self._model, "ecutwfc"),
            lambda cutoffs: max(cutoffs[0]),
        )
        ipw.dlink(
            (self._model, "cutoffs"),
            (self._model, "ecutrho"),
            lambda cutoffs: max(cutoffs[1]),
        )

    def render(self):
        if self.rendered:
            return

        self.family_prompt = ipw.HTML()

        self.family_help = ipw.HTML()
        ipw.dlink(
            (self._model, "family_help_message"),
            (self.family_help, "value"),
        )

        self.functional_prompt = ipw.HTML("""
            <div class="pseudo-text"><b>
                Exchange-correlation  functional</b>
            </div>
        """)

        self.functional_help = ipw.HTML("""
            <div class="pseudo-text">
                The exchange-correlation energy is calculated using this functional. We
                currently provide support for two well-established generalized gradient
                approximation (GGA) functionals: PBE and PBEsol.
            </div>
        """)

        self.functional = ipw.Dropdown(style={"description_width": "initial"})
        ipw.dlink(
            (self._model, "functional_options"),
            (self.functional, "options"),
        )
        ipw.link(
            (self._model, "functional"),
            (self.functional, "value"),
        )
        ipw.dlink(
            (self._model, "override"),
            (self.functional, "disabled"),
            lambda override: not override,
        )

        self.library = ipw.ToggleButtons(layout=ipw.Layout(max_width="80%"))
        ipw.dlink(
            (self._model, "library_options"),
            (self.library, "options"),
        )
        ipw.link(
            (self._model, "library"),
            (self.library, "value"),
        )
        ipw.dlink(
            (self._model, "override"),
            (self.library, "disabled"),
            lambda override: not override,
        )

        self.setter_widget_helper = ipw.HTML("""
            <div class="pseudo-text">
                The pseudopotential for each kind of atom in the structure can be
                custom set. The default pseudopotential and cutoffs are get from
                the pseudo family. The cutoffs used for the calculation are the
                maximum of the default from all pseudopotentials and can be custom
                set.
            </div>
        """)

        self.setter_widget = ipw.VBox()

        self._status_message = StatusHTML(clear_after=20)
        ipw.dlink(
            (self._model, "status_message"),
            (self._status_message, "message"),
        )

        self.cutoff_helper = ipw.HTML("""
            <div class="pseudo-text">
                Please set the cutoffs for the calculation. The default cutoffs are get
                from the pseudo family.
            </div>
        """)
        self.ecutwfc = ipw.FloatText(
            description="Wavefunction cutoff (Ry)",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "ecutwfc"),
            (self.ecutwfc, "value"),
        )
        ipw.dlink(
            (self._model, "override"),
            (self.ecutwfc, "disabled"),
            lambda override: not override,
        )
        self.ecutrho = ipw.FloatText(
            description="Charge density cutoff (Ry)",
            style={"description_width": "initial"},
        )
        ipw.link(
            (self._model, "ecutrho"),
            (self.ecutrho, "value"),
        )
        ipw.dlink(
            (self._model, "override"),
            (self.ecutrho, "disabled"),
            lambda override: not override,
        )

        self.container = ipw.VBox()

        self.children = [
            ipw.HTML("""
                <div style="padding-top: 0px; padding-bottom: 10px">
                    <h4>Accuracy and precision</h4>
                </div>
            """),
            ipw.HBox(
                children=[
                    ipw.HTML(
                        """
                        <div class="pseudo-text">
                            The exchange-correlation functional and pseudopotential
                            library is set by the <b>protocol</b> configured in the
                            "Workflow" tab. Here you can override the defaults if
                            desired.
                        </div>
                        """,
                        layout=ipw.Layout(max_width="60%"),
                    ),
                ],
                layout=ipw.Layout(height="50px", justify_content="space-between"),
            ),
            ipw.HBox(
                [
                    ipw.VBox(
                        children=[
                            self.functional_prompt,
                            self.functional,
                            self.functional_help,
                        ],
                        layout=ipw.Layout(max_width="40%"),
                    ),
                    ipw.VBox(
                        children=[
                            self.family_prompt,
                            self.library,
                            self.family_help,
                        ],
                        layout=ipw.Layout(max_width="60%"),
                    ),
                ]
            ),
            self.container,
            self.cutoff_helper,
            ipw.HBox(
                children=[
                    self.ecutwfc,
                    self.ecutrho,
                ],
            ),
            self._status_message,
        ]

        self.rendered = True

        self.refresh()

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_protocol_change(self, _):
        self.refresh(specific="protocol")

    def _on_spin_orbit_change(self, _):
        self._model.update_library_options()

    def _on_override_change(self, change):
        super()._on_override_change(change)
        self._toggle_setter_widgets()

    def _on_family_parameters_change(self, _):
        self._model.update_family()

    def _on_family_change(self, _):
        self._update_family_link()
        self._model.update_default_pseudos()
        self._model.update_default_cutoffs()

    def _update(self, specific=""):
        if self.updated:
            return
        self._show_loading()
        if not self._model.loaded_from_process:
            self._model.update(specific)
        self._build_setter_widgets()
        self._toggle_setter_widgets()
        self._model.update_library_options()
        self._update_family_link()
        self.updated = True

    def _update_family_link(self):
        if not self.rendered:
            return

        library, accuracy = self._model.library.split()
        if library == "SSSP":
            pseudo_family_link = (
                f"https://www.materialscloud.org/discover/sssp/table/{accuracy}"
            )
        else:
            pseudo_family_link = "http://www.pseudo-dojo.org/"

        self.family_prompt.value = f"""
            <div class="pseudo-text">
                <b>
                    <a href="{pseudo_family_link}" target="_blank">
                        Pseudopotential family
                    </a>
                </b>
            </div>
        """

    def _show_loading(self):
        if self.rendered:
            self.setter_widget.children = [self.loading_message]

    def _build_setter_widgets(self):
        if not self.rendered:
            return

        children = []

        kinds = self._model.input_structure.kinds if self._model.input_structure else []

        for index, kind in enumerate(kinds):
            upload_widget = PseudoUploadWidget(kind_name=kind.name)
            pseudo_link = ipw.link(
                (self._model, "dictionary"),
                (upload_widget, "pseudo"),
                [
                    lambda dictionary, kind_name=kind.name: orm.load_node(
                        dictionary.get(kind_name)
                    ),
                    lambda pseudo, kind_name=kind.name: {
                        **self._model.dictionary,
                        kind_name: pseudo.uuid,
                    },
                ],
            )
            cutoffs_link = ipw.dlink(
                (self._model, "cutoffs"),
                (upload_widget, "cutoffs"),
                lambda cutoffs, index=index: [cutoffs[0][index], cutoffs[1][index]]
                if len(cutoffs[0]) > index
                else [0.0, 0.0],
            )
            upload_widget.render()

            self.links.extend(
                [
                    pseudo_link,
                    cutoffs_link,
                    *upload_widget.links,
                ]
            )

            children.append(upload_widget)

        self.setter_widget.children = children

    def _toggle_setter_widgets(self):
        if not self.rendered:
            return
        if self._model.override:
            self.container.children = [
                self.setter_widget_helper,
                self.setter_widget,
            ]
        else:
            self.container.children = []


# TODO implement/improve MVC in this widget
class PseudoUploadWidget(ipw.HBox):
    """Class that allows to upload pseudopotential from user's computer."""

    pseudo = tl.Instance(UpfData, allow_none=True)
    cutoffs = tl.List(tl.Float(), [])
    error_message = tl.Unicode(allow_none=True)

    def __init__(self, kind_name, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading pseudopotential uploader")],
            **kwargs,
        )

        self.kind_name = kind_name

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.file_upload = ipw.FileUpload(
            description="Upload",
            multiple=False,
        )
        self.pseudo_text = ipw.Text(description=self.kind_name)
        self.file_upload.observe(self._on_file_upload, "value")

        cutoffs_message_template = """
            <div class="pseudo-text">
                Recommended ecutwfc: <b>{ecutwfc} Ry</b> ecutrho: <b>{ecutrho} Ry</b>
            </div>
        """

        self.cutoff_message = ipw.HTML()

        pseudo_link = ipw.dlink(
            (self, "pseudo"),
            (self.pseudo_text, "value"),
            lambda pseudo: pseudo.filename if pseudo else "",
        )

        cutoff_link = ipw.dlink(
            (self, "cutoffs"),
            (self.cutoff_message, "value"),
            lambda cutoffs: cutoffs_message_template.format(
                ecutwfc=cutoffs[0] if len(cutoffs) else "not set",
                ecutrho=cutoffs[1] if len(cutoffs) else "not set",
            ),
        )

        self.links = [pseudo_link, cutoff_link]

        self.error_message = None

        self.children = [
            self.pseudo_text,
            self.file_upload,
            self.cutoff_message,
        ]

        self.rendered = True

    def _on_file_upload(self, change=None):
        """When file upload button is pressed."""
        filename, item = next(iter(change["new"].items()))
        content = item["content"]

        # Order matters make sure when pseudo change
        # the pseudo_filename is set
        with self.hold_trait_notifications():
            self.pseudo = UpfData(io.BytesIO(content), filename=filename)
            self.pseudo.store()

            # check if element is matched with the pseudo
            element = "".join([i for i in self.kind_name if not i.isdigit()])
            if element != self.pseudo.element:
                self.error_message = f"""<div class="alert alert-danger"> ERROR: Element {self.kind_name} is not matched with the pseudo {self.pseudo.element}</div>"""
                self._reset()
            else:
                self.pseudo_text.value = filename

    def _reset(self):
        """Reset the widget to the initial state."""
        self.pseudo = None
        self.cutoffs = []
