from __future__ import annotations

import io
from copy import deepcopy

import ipywidgets as ipw
import traitlets as tl
from aiida_pseudo.common.units import U
from upf_tools import UPFDict

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


_PSEUDO_ALERT_TEMPLATE = """
    <div class="alert alert-{alert_type} pseudo-warning">
        {message}
    </div>
"""


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
            value=_PSEUDO_ALERT_TEMPLATE.format(
                alert_type="warning",
                message="""
                        You have uploaded a custom pseudopotential.
                        <br>
                        <b>Please make sure to use the same functional for all pseudopotentials.</b>
                        <br>
                        To reset to the protocol-derived pseudopotentials, please select a functional and a family.
                    """,
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
                    self.family_prompt,
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
        self._update_family_link()
        if not (self._model.library and self._model.functional):
            return
        self._model.show_upload_warning = False
        self._model.update_family()

    def _on_family_change(self, _):
        if not self._model.family:
            return
        self._model.update_default_pseudos()
        self._model.update_default_cutoffs()

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
        self._update_family_link()
        self.updated = True

    def _update_family_link(self):
        if not self.rendered:
            return

        # Remove link if library is unset
        if not self._model.library:
            self.family_prompt.value = "<h4>Pseudopotential family</h4>"
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
            upload_widget = PseudoUploadWidget(kind.name, kind.symbol)

            def on_default_pseudo(
                _=None,
                kind_name=kind.name,
                widget=upload_widget,
            ):
                try:
                    pseudo = orm.load_node(self._model.dictionary.get(kind_name))
                except Exception:
                    pseudo = None
                widget.pseudo = pseudo
                widget.uploaded = False

            self._model.observe(
                on_default_pseudo,
                "dictionary",
            )

            def on_pseudo_upload(
                change,
                kind_name=kind.name,
                index=index,
                widget=upload_widget,
            ):
                if not change["new"]:
                    return

                self._model.functional = None
                self._model.library = None
                self._model.family = None
                self._model.show_upload_warning = True

                self._model.dictionary = {
                    **self._model.dictionary,
                    kind_name: widget.pseudo.uuid,
                }

                cutoffs: list = deepcopy(self._model.cutoffs)  # type: ignore
                cutoffs[0][index] = widget.cutoffs[0]
                cutoffs[1][index] = widget.cutoffs[1]
                self._model.cutoffs = cutoffs

            upload_widget.observe(
                on_pseudo_upload,
                "uploaded",
            )

            cutoffs_link = ipw.dlink(
                (self._model, "cutoffs"),
                (upload_widget, "cutoffs"),
                lambda cutoffs, index=index: [cutoffs[0][index], cutoffs[1][index]]
                if len(cutoffs[0]) > index
                else [0.0, 0.0],
            )

            upload_widget.render()

            on_default_pseudo()

            self.links.extend(
                [
                    cutoffs_link,
                    *upload_widget.links,
                ]
            )

            children.append(upload_widget)

        self.setter_widget.children = children


class PseudoUploadWidget(ipw.VBox):
    """Class that allows to upload pseudopotential from user's computer."""

    pseudo = tl.Instance(UpfData, allow_none=True)
    cutoffs = tl.List(tl.Float(), [])
    uploaded = tl.Bool(False)

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

        self.message_box = StatusHTML(clear_after=5)

        self.pseudo_row = ipw.HBox(
            children=[
                self.pseudo_text,
                self.file_upload,
                self.cutoff_message,
            ]
        )

        self.children = [
            self.pseudo_row,
            self.message_box,
        ]

        self.rendered = True

    def _on_file_upload(self, change=None):
        if not change or not change["new"]:
            return

        filename, item = next(iter(change["new"].items()))
        content = item["content"]

        try:
            uploaded_pseudo = UpfData(io.BytesIO(content), filename=filename)
        except Exception:
            self.message_box.message = _PSEUDO_ALERT_TEMPLATE.format(
                alert_type="danger",
                message=f"{filename} is not a valid UPF file",
            )
            return

        # Wrong element
        if uploaded_pseudo.element != self.kind_symbol:
            self.message_box.message = _PSEUDO_ALERT_TEMPLATE.format(
                alert_type="danger",
                message=f"Pseudo element {uploaded_pseudo.element} does not match {self.kind_symbol}",
            )
            self._reset_uploader()
            return

        # Existing pseudo
        if existing_pseudo := (
            orm.QueryBuilder()
            .append(
                UpfData,
                filters={"attributes.md5": uploaded_pseudo.md5},
            )
            .first(flat=True)
        ):
            uploaded_pseudo = existing_pseudo
            message = _PSEUDO_ALERT_TEMPLATE.format(
                alert_type="info",
                message=f"Identical pseudo detected. Loading pseudo (UUID={uploaded_pseudo.uuid})",
            )

        # New pseudo but existing filename
        elif (
            orm.QueryBuilder()
            .append(
                UpfData,
                filters={"attributes.filename": filename},
            )
            .count()
        ):
            self.message_box.message = _PSEUDO_ALERT_TEMPLATE.format(
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
            message = _PSEUDO_ALERT_TEMPLATE.format(
                alert_type="success",
                message=f"{filename} uploaded successfully",
            )

        if pseudo_family := self._get_pseudo_family(uploaded_pseudo.md5):
            current_unit = pseudo_family.get_cutoffs_unit()
            cutoff_dict = pseudo_family.get_cutoffs()
            cutoff = {
                key: U.Quantity(v, current_unit).to("Ry").to_tuple()[0]
                for key, v in cutoff_dict.get(self.kind_symbol, {}).items()
            }
            cutoffs = [
                cutoff.get("cutoff_wfc", 0.0),
                cutoff.get("cutoff_rho", 0.0),
            ]
        else:
            with uploaded_pseudo.as_path() as pseudo_path:
                upf_dict = UPFDict.from_upf(pseudo_path.as_posix())
            cutoffs = [
                float(int(upf_dict["header"].get("wfc_cutoff", 0))),
                float(int(upf_dict["header"].get("rho_cutoff", 0))),
            ]

        self.cutoffs = cutoffs

        self.message_box.message = message
        self.pseudo = uploaded_pseudo
        self.uploaded = True
        self._reset_uploader()

    def _get_pseudo_family(self, pseudo_md5):
        return (
            orm.QueryBuilder()
            .append(
                UpfData,
                filters={"attributes.md5": pseudo_md5},
                tag="pseudo",
            )
            .append(
                (
                    SsspFamily,
                    PseudoDojoFamily,
                ),
                with_node="pseudo",
            )
            .first(flat=True)
        )

    def _reset_uploader(self):
        self.file_upload.metadata = []
        self.file_upload.data = []
        self.file_upload._counter = 0
