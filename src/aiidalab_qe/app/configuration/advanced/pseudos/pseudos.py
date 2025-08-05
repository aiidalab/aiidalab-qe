from __future__ import annotations

import io
from copy import deepcopy

import ipywidgets as ipw
import traitlets as tl

from aiida import orm
from aiidalab_qe.common.widgets import HBoxWithUnits, LoadingWidget
from aiidalab_qe.utils import (
    UpfData,
    generate_alert,
    get_pseudo_by_filename,
    get_pseudo_by_md5,
    populate_extras,
)
from aiidalab_widgets_base.utils import StatusHTML

from ..subsettings import AdvancedConfigurationSubSettingsPanel
from .model import PseudosConfigurationSettingsModel


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
            self._on_functionals_change,
            "functionals",
        )
        self._model.observe(
            self._on_family_change,
            "family",
        )
        self._model.observe(
            self._on_dictionary_change,
            "dictionary",
        )
        self._model.observe(
            self._on_cutoffs_change,
            "cutoffs",
        )
        self._model.observe(
            self._on_show_upload_warning_change,
            "show_upload_warning",
        )

    def render(self):
        if self.rendered:
            return

        self.family_header = ipw.HTML()
        ipw.dlink(
            (self._model, "family_header"),
            (self.family_header, "value"),
        )

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

        self._status_message = StatusHTML(clear_after=20)
        ipw.link(
            (self._model, "status_message"),
            (self._status_message, "message"),
        )

        self._warning_message = ipw.HTML(
            value=generate_alert(
                alert_type="warning",
                message="""
                    <h4><b>⚠️ Pseudopotential upload detected ⚠️</b></h4>
                    Adjust plane-wave cutoff energies below to ensure convergence
                    <br>
                    To reset to protocol-derived pseudopotentials, select both an
                    exchange-correlation functional and a pseudopotential family above
                """,
                class_="text-center",
                style_="line-height: 2;",
            ),
            layout=ipw.Layout(
                display="block" if self._model.show_upload_warning else "none",
            ),
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
                    self.family_header,
                    self.library,
                    self.family_help,
                ],
            ),
            ipw.HTML("<h4>Pseudopotentials</h4>"),
            self._warning_message,
            ipw.HTML("""
                <div class="pseudo-text">
                    The pseudopotential for each kind of atom in the structure can be
                    custom set.
                    <br>
                    The default pseudopotential and cutoffs are taken from the
                    pseudopotential family.
                    <br>
                    Pseudopotential information is shown (if provided) in the following order:
                    <ul>
                        <li>
                            Recommended wavefunction cutoff energy (Ry)
                        </li>
                        <li>
                            Recommended charge density cutoff energy (Ry)
                        </li>
                        <li>
                            Exchange-correlation functional
                        </li>
                        <li>
                            Relativistic treatment
                        </li>
                    </ul>
                </div>
            """),
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
        self._model.update_family_header()
        self._model.update_functionals()
        self._model.update_blockers()

    def _on_functionals_change(self, _):
        self._model.functional = (
            self._model.functionals[0]
            if len(set(self._model.functionals)) == 1
            else None
        )
        self._model.update_blockers()

    def _on_family_change(self, _):
        self._model.show_upload_warning = not self._model.family
        self._model.update_dictionary()

    def _on_dictionary_change(self, _):
        self._model.update_cutoffs()

    def _on_cutoffs_change(self, change):
        cutoffs = change["new"]  # [[ecutwfc...], [ecutrho...]]
        self._model.ecutwfc = max(cutoffs[0])
        self._model.ecutrho = max(cutoffs[1])

    def _on_show_upload_warning_change(self, change):
        if not self.rendered:
            return
        self._warning_message.layout.display = "block" if change["new"] else "none"

    def _update(self, specific=""):
        if self.updated:
            return
        self._show_loading()
        if not self._model.loaded_from_process or (specific and specific != "widgets"):
            self._model.update(specific)
        self._build_setter_widgets()
        self._model.update_library_options()
        self._model.update_family_header()
        self.updated = True

    def _show_loading(self):
        if self.rendered:
            self.setter_widget.children = [self.loading_message]

    def _build_setter_widgets(self):
        if not self.rendered:
            return

        children = []

        kinds = self._model.input_structure.kinds if self._model.has_structure else []

        for index, kind in enumerate(kinds):
            upload_widget = PseudoUploadWidget(kind.name, kind.symbol)

            def on_default_pseudo(
                _=None,
                widget: PseudoUploadWidget = upload_widget,
                kind_name=kind.name,
                index=index,
            ):
                if not self._model.dictionary:
                    return
                try:
                    pseudo = orm.load_node(self._model.dictionary.get(kind_name))
                except Exception:
                    pseudo = None
                widget.pseudo = pseudo
                widget.cutoffs = (
                    [
                        self._model.cutoffs[0][index],
                        self._model.cutoffs[1][index],
                    ]
                    if len(self._model.cutoffs[0]) > index
                    else [0.0, 0.0]
                )
                widget.update_pseudo_info()
                widget.uploaded = False

            self._model.observe(
                on_default_pseudo,
                "dictionary",
            )

            def on_pseudo_upload(
                change,
                widget: PseudoUploadWidget = upload_widget,
                kind_name=kind.name,
                index=index,
            ):
                if not (change["new"] and widget.pseudo):
                    return

                self._model.family = None
                self._model.library = None
                self._model.functional = None

                self._model.dictionary = {
                    **self._model.dictionary,
                    kind_name: widget.pseudo.uuid,
                }

                functional = widget.pseudo.base.extras.get("functional", None)
                functionals = [*self._model.functionals]
                functionals[index] = functional
                # The following double-setting is done to force the blockers check,
                # which now also bail early if the functionals are empty.
                self._model.functionals = []
                self._model.functionals = functionals

                cutoffs: list = deepcopy(self._model.cutoffs)  # type: ignore
                cutoffs[0][index] = widget.cutoffs[0]
                cutoffs[1][index] = widget.cutoffs[1]
                self._model.cutoffs = cutoffs

            upload_widget.observe(
                on_pseudo_upload,
                "uploaded",
            )

            upload_widget.render()

            on_default_pseudo()

            self.links.extend(upload_widget.links)

            children.append(upload_widget)

        self.setter_widget.children = children


# TODO refactor as MVC
class PseudoUploadWidget(ipw.VBox):
    """Class that allows to upload pseudopotential from user's computer."""

    pseudo = tl.Instance(UpfData, allow_none=True)
    filename = tl.Unicode(allow_none=True)
    info = tl.Unicode(allow_none=True)
    cutoffs = tl.List(tl.Float(), [])
    uploaded = tl.Bool(False)
    message = tl.Unicode(allow_none=True)

    PSEUDO_INFO_MESSAGE = """
        <div class="pseudo-text">
            {ecutwfc} | {ecutrho} | {functional} | {relativistic}
        </div>
    """

    def __init__(self, kind_name, kind_symbol, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading pseudopotential uploader")],
            **kwargs,
        )

        self.kind_name = kind_name
        self.kind_symbol = kind_symbol

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.pseudo_filename = ipw.Text(
            description=self.kind_name,
            style={"description_width": "50px"},
            disabled=True,
        )
        filename_link = ipw.dlink(
            (self, "pseudo"),
            (self.pseudo_filename, "value"),
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

        self.pseudo_info = ipw.HTML()
        info_link = ipw.dlink(
            (self, "info"),
            (self.pseudo_info, "value"),
        )

        self.message_box = StatusHTML(clear_after=5)
        message_link = ipw.link(
            (self, "message"),
            (self.message_box, "message"),
        )

        self.links = [
            filename_link,
            info_link,
            message_link,
        ]

        self.pseudo_row = ipw.HBox(
            children=[
                self.pseudo_filename,
                self.file_upload,
                self.pseudo_info,
            ]
        )

        self.children = [
            self.pseudo_row,
            self.message_box,
        ]

        self.rendered = True

    def update_pseudo_info(self):
        if not self.pseudo:
            self.info = ""
            return
        self.info = self.PSEUDO_INFO_MESSAGE.format(
            ecutwfc=self.cutoffs[0] or "?",
            ecutrho=self.cutoffs[1] or "?",
            functional=self.pseudo.base.extras.get("functional", None) or "?",
            relativistic=self.pseudo.base.extras.get("relativistic", None) or "N/A",
        )

    @tl.observe("pseudo")
    def _on_pseudo_change(self, _):
        if not self.pseudo:
            return
        populate_extras(self.pseudo)

    def _on_file_upload(self, change=None):
        if not change or not change["new"]:
            return

        filename, item = next(iter(change["new"].items()))
        content = item["content"]

        try:
            uploaded_pseudo = UpfData(io.BytesIO(content), filename=filename)
        except Exception:
            self.message = generate_alert(
                alert_type="danger",
                message=f"{filename} is not a valid UPF file",
            )
            self._reset_uploader()
            return

        # Wrong element
        if uploaded_pseudo.element != self.kind_symbol:
            self.message = generate_alert(
                alert_type="danger",
                message=f"Pseudo element {uploaded_pseudo.element} does not match {self.kind_symbol}",
            )
            self._reset_uploader()
            return

        # Existing pseudo
        if existing_pseudo := get_pseudo_by_md5(uploaded_pseudo.md5):
            uploaded_pseudo = existing_pseudo
            message = generate_alert(
                alert_type="info",
                message=f"Identical pseudo detected. Loading pseudo (UUID={uploaded_pseudo.uuid})",
            )

        # New pseudo but existing filename
        elif get_pseudo_by_filename(filename):
            self.message = generate_alert(
                alert_type="warning",
                message=f"""
                    {filename} found in database with different content.
                    <br>
                    Please rename your file before uploading.
                """,
            )
            self._reset_uploader()
            return

        # Valid new pseudo
        else:
            uploaded_pseudo.store()
            message = generate_alert(
                alert_type="success",
                message=f"{filename} uploaded successfully",
            )

        self.pseudo = uploaded_pseudo
        self.message = message
        self.uploaded = True
        self._reset_uploader()

    def _reset_uploader(self):
        self.file_upload.metadata = []
        self.file_upload.data = []
        self.file_upload._counter = 0
