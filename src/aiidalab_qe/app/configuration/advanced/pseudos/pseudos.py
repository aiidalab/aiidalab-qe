from __future__ import annotations

import io

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiida.plugins import DataFactory, GroupFactory
from aiidalab_qe.common.widgets import HBoxWithUnits, LoadingWidget
from aiidalab_widgets_base.utils import StatusHTML

from ..subsettings import AdvancedConfigurationSubSettingsPanel
from .model import PseudosConfigurationSettingsModel

UpfData = DataFactory("pseudo.upf")
SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")


class PseudosConfigurationSettingsPanel(
    AdvancedConfigurationSubSettingsPanel[PseudosConfigurationSettingsModel],
):
    def __init__(self, model: PseudosConfigurationSettingsModel, **kwargs):
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
        self._model.observe(
            self._on_reset_trigger_change,
            "pseudo_filename_reset_trigger",
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

        self.functional_help = ipw.HTML("""
            <div class="pseudo-text">
                The exchange-correlation energy is calculated using this functional.
                <br>
                We currently provide support for two well-established generalized
                gradient approximation (GGA) functionals: PBE and PBEsol.
            </div>
        """)

        self.functional = ipw.ToggleButtons()
        ipw.dlink(
            (self._model, "functional_options"),
            (self.functional, "options"),
        )
        ipw.link(
            (self._model, "functional"),
            (self.functional, "value"),
        )

        self.library = ipw.ToggleButtons()
        ipw.dlink(
            (self._model, "library_options"),
            (self.library, "options"),
        )
        ipw.link(
            (self._model, "library"),
            (self.library, "value"),
        )

        self.setter_widget = ipw.VBox()

        self._status_message = StatusHTML(clear_after=20)
        ipw.dlink(
            (self._model, "status_message"),
            (self._status_message, "message"),
        )

        self.ecutwfc = ipw.FloatText(
            description="Wavefunction",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "ecutwfc"),
            (self.ecutwfc, "value"),
        )
        self.ecutrho = ipw.FloatText(
            description="Charge density",
            style={"description_width": "150px"},
        )
        ipw.link(
            (self._model, "ecutrho"),
            (self.ecutrho, "value"),
        )

        self.children = [
            ipw.HTML("<h2>Accuracy and precision</h2>"),
            ipw.HTML("""
                <div class="pseudo-text">
                    The exchange-correlation functional and pseudopotential library is
                    set by the <b>protocol</b> configured in the <b>Basic settings</b>
                    tab.
                    <br>
                    Here you can override the defaults if desired.
                </div>
            """),
            ipw.VBox(
                children=[
                    ipw.HTML("<h4>Exchange-correlation functional</h4>"),
                    self.functional,
                    self.functional_help,
                ],
            ),
            ipw.VBox(
                children=[
                    self.family_prompt,
                    self.library,
                    self.family_help,
                ],
            ),
            ipw.HTML("<h4>Pseudopotentials</h4>"),
            ipw.HTML("""
                <div class="pseudo-text">
                    The pseudopotential for each kind of atom in the structure can be
                    custom set.
                    <br>
                    The default pseudopotential and cutoffs are taken from the
                    pseudopotential family.
                    <br>
                    Recommended wavefunction (ψ) and charge density (ρ) cutoffs are
                    given to the right of each pseudopotential.
                </div>
            """),  # noqa: RUF001
            self.setter_widget,
            ipw.HTML("<h4>Cutoffs</h4>"),
            ipw.HTML("""
                <div style="line-height: 1.4;">
                    The
                    <a
                        href="https://www.quantum-espresso.org/Doc/INPUT_PW.html#idm312"
                        target="_blank"
                    >
                        default cutoffs
                    </a> used for the calculation are the maximum of the default cutoffs
                    from all pseudopotentials.
                    <br>
                    You can override them here.
                </div>
            """),
            HBoxWithUnits(self.ecutwfc, "Ry"),
            HBoxWithUnits(self.ecutrho, "Ry"),
            self._status_message,
        ]

        self.rendered = True

        self.refresh(specific="widgets")

    def _on_input_structure_change(self, _):
        self.refresh(specific="structure")

    def _on_protocol_change(self, _):
        self.refresh(specific="protocol")

    def _on_spin_orbit_change(self, _):
        self._model.update_library_options()

    def _on_family_parameters_change(self, _):
        self._model.update_family()

    def _on_family_change(self, _):
        self._update_family_link()
        self._model.update_default_pseudos()
        self._model.update_default_cutoffs()

    def _on_reset_trigger_change(self, _):
        if not self.rendered:
            return
        upload_widget: PseudoUploadWidget
        for upload_widget in self.setter_widget.children:
            upload_widget.reset_pseudo_filename_widget()

    def _update(self, specific=""):
        if self.updated:
            return
        self._show_loading()
        if not self._model.loaded_from_process or (specific and specific != "widgets"):
            self._model.update(specific)
        self._build_setter_widgets()
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
            <h4>
                <a href="{pseudo_family_link}" target="_blank">
                    Pseudopotential family
                </a>
            </h4>
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

        self.pseudo_text = ipw.Text(
            description=self.kind_name,
            style={"description_width": "50px"},
        )
        pseudo_link = ipw.dlink(
            (self, "pseudo"),
            (self.pseudo_text, "value"),
            lambda pseudo: pseudo.filename if pseudo else "",
        )
        self.file_upload = ipw.FileUpload(
            description="Upload",
            multiple=False,
            layout=ipw.Layout(
                width="fit-content",
                margin="2px 10px 2px 2px",
            ),
        )
        self.file_upload.observe(self._on_file_upload, "value")

        cutoffs_message_template = """
            <div class="pseudo-text">
                ψ: <b>{ecutwfc} Ry</b> | ρ: <b>{ecutrho} Ry</b>
            </div>
        """  # noqa: RUF001

        self.cutoff_message = ipw.HTML()
        cutoff_link = ipw.dlink(
            (self, "cutoffs"),
            (self.cutoff_message, "value"),
            lambda cutoffs: cutoffs_message_template.format(
                ecutwfc=cutoffs[0] if len(cutoffs) else "not set",
                ecutrho=cutoffs[1] if len(cutoffs) else "not set",
            ),
        )

        self.links = [
            pseudo_link,
            cutoff_link,
        ]

        self.error_message = None

        self.children = [
            self.pseudo_text,
            self.file_upload,
            self.cutoff_message,
        ]

        self.rendered = True

    def reset_pseudo_filename_widget(self):
        self.pseudo_text.value = self.pseudo.filename if self.pseudo else ""

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
