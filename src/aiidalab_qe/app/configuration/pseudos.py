from __future__ import annotations

import io

import ipywidgets as ipw
import traitlets as tl
from aiida_pseudo.common.units import U

from aiida import orm
from aiida.common import exceptions
from aiida.plugins import DataFactory, GroupFactory
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_qe.common.widgets import LoadingWidget
from aiidalab_qe.setup.pseudos import (
    PSEUDODOJO_VERSION,
    SSSP_VERSION,
    PseudoFamily,
)
from aiidalab_widgets_base.utils import StatusHTML

from .model import config_model as model

UpfData = DataFactory("pseudo.upf")
SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")


class PseudoSettings(ipw.VBox):
    """Widget to set the pseudopotentials for the calculation."""

    protocol = tl.Unicode(allow_none=True)
    input_structure = tl.Instance(klass=orm.StructureData, allow_none=True)
    spin_orbit = tl.Unicode()
    pseudo_family = tl.Unicode()

    def __init__(self, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading pseudopotentials widget")],
            **kwargs,
        )

        self.setter_widget_links = []

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.PSEUDO_HELP_SOC = """
            <div class="pseudo-text">
                Spin-orbit coupling (SOC) calculations are supported exclusively with
                PseudoDojo pseudopotentials. PseudoDojo offers these pseudopotentials
                in two versions: standard and stringent. Here, we utilize the FR
                (fully relativistic) type from PseudoDojo. Please ensure you choose
                appropriate cutoff values for your calculations.
            </div>
        """

        self.PSEUDO_HELP_WO_SOC = """
            <div class="pseudo-text">
                If you are unsure, select 'SSSP efficiency', which for most
                calculations will produce sufficiently accurate results at
                comparatively small computational costs. If your calculations require a
                higher accuracy, select 'SSSP accuracy' or 'PseudoDojo stringent',
                which will be computationally more expensive. SSSP is the standard
                solid-state pseudopotentials. The PseudoDojo used here has the SR
                relativistic type.
            </div>
        """

        # XXX: the link is not correct after add pseudo dojo
        self.family_prompt = ipw.HTML()

        self.family_help = ipw.HTML(self.PSEUDO_HELP_WO_SOC)

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

        self.override = ipw.Checkbox(
            description="",
            indent=False,
            layout=ipw.Layout(max_width="10%"),
        )
        ipw.link(
            (model.pseudos, "override"),
            (self.override, "value"),
        )
        ipw.dlink(
            (model, "input_structure"),
            (self.override, "disabled"),
            lambda structure: structure is None,
        )
        self.override.observe(self._on_override_change, "value")

        self.functional = ipw.Dropdown(
            options=["PBE", "PBEsol"],
            style={"description_width": "initial"},
        )
        ipw.link(
            (model.pseudos, "functional"),
            (self.functional, "value"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.functional, "disabled"),
            lambda override: not override,
        )
        self.functional.observe(self._on_family_parameters_change, "value")

        self.library = ipw.ToggleButtons(
            options=[
                "SSSP efficiency",
                "SSSP precision",
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ],
            layout=ipw.Layout(max_width="80%"),
        )
        ipw.link(
            (model.pseudos, "library"),
            (self.library, "value"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.library, "disabled"),
            lambda override: not override,
        )
        self.library.observe(self._on_family_parameters_change, "value")

        self.setter_widget = ipw.VBox()

        self._status_message = StatusHTML(clear_after=20)

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
            (model.pseudos, "ecutwfc"),
            (self.ecutwfc, "value"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.ecutwfc, "disabled"),
            lambda override: not override,
        )
        self.ecutrho = ipw.FloatText(
            description="Charge density cutoff (Ry)",
            style={"description_width": "initial"},
        )
        ipw.link(
            (model.pseudos, "ecutrho"),
            (self.ecutrho, "value"),
        )
        ipw.dlink(
            (self.override, "value"),
            (self.ecutrho, "disabled"),
            lambda override: not override,
        )

        self.container = ipw.VBox(
            children=[
                ipw.HTML("""
                <div class="pseudo-text">
                    The pseudopotential for each kind of atom in the structure can be
                    custom set. The default pseudopotential and cutoffs are get from
                    the pseudo family. The cutoffs used for the calculation are the
                    maximum of the default from all pseudopotentials and can be custom
                    set.
                </div>
                """),
            ]
        )

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
                    ipw.HBox(
                        children=[
                            ipw.HTML(
                                value="<b>Override</b>",
                                layout=ipw.Layout(margin="0 5px 0 0"),
                            ),
                            self.override,
                        ],
                        layout=ipw.Layout(max_width="20%"),
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

        ipw.dlink(
            (model.pseudos, "family"),
            (self, "pseudo_family"),
        )
        ipw.dlink(
            (model, "protocol"),
            (self, "protocol"),
        )
        ipw.dlink(
            (model, "input_structure"),
            (self, "input_structure"),
        )
        ipw.dlink(
            (model, "spin_orbit"),
            (self, "spin_orbit"),
        )

        self.rendered = True

    def reset(self):
        self._unsubscribe_setter_links()
        model.pseudos.reset()

    @tl.observe("protocol")
    def _on_protocol_change(self, _=None):
        self._update_family_parameters()

    @tl.observe("input_structure")
    def _on_input_structure_change(self, _=None):
        self._build_setter_widgets()

    @tl.observe("spin_orbit")
    def _on_spin_orbit_change(self, _):
        self._update_library_options()

    @tl.observe("pseudo_family")
    def _on_pseudo_family_change(self, _=None):
        self._update_family_link()
        self._update_pseudos()
        self._update_cutoffs()

    def _on_override_change(self, change):
        self._toggle_setter_widgets(change)
        if not change["new"]:
            self.reset()

    def _on_family_parameters_change(self, _=None):
        self._update_family()

    def _update_family(self):
        """Set pseudo family according to the selected library and functional."""
        library, accuracy = model.pseudos.library.split()
        functional = model.pseudos.functional
        # XXX (jusong.yu): a validator is needed to check the family string is consistent with the list of pseudo families defined in the setup_pseudos.py
        if library == "PseudoDojo":
            if model.spin_orbit == "soc":
                pseudo_family_string = (
                    f"PseudoDojo/{PSEUDODOJO_VERSION}/{functional}/FR/{accuracy}/upf"
                )
            else:
                pseudo_family_string = (
                    f"PseudoDojo/{PSEUDODOJO_VERSION}/{functional}/SR/{accuracy}/upf"
                )
        elif library == "SSSP":
            pseudo_family_string = f"SSSP/{SSSP_VERSION}/{functional}/{accuracy}"
        else:
            raise ValueError(
                f"Unknown pseudo family parameters: {library} | {accuracy}"
            )

        model.pseudos.family = pseudo_family_string

    def _update_family_parameters(self):
        if model.spin_orbit == "soc":
            if model.protocol in ["fast", "moderate"]:
                pseudo_family_string = "PseudoDojo/0.4/PBE/FR/standard/upf"
            else:
                pseudo_family_string = "PseudoDojo/0.4/PBE/FR/stringent/upf"
        else:
            pseudo_family_string = PwBaseWorkChain.get_protocol_inputs(model.protocol)[
                "pseudo_family"
            ]

        pseudo_family = PseudoFamily.from_string(pseudo_family_string)

        with self.hold_trait_notifications():
            model.pseudos.library = f"{pseudo_family.library} {pseudo_family.accuracy}"
            model.pseudos.functional = pseudo_family.functional

    def _update_library_options(self):
        """Update pseudo library selection options w.r.t spin orbit."""
        if model.spin_orbit == "soc":
            self.library.options = [
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ]
            self.family_help.value = self.PSEUDO_HELP_SOC
        else:
            self.library.options = [
                "SSSP efficiency",
                "SSSP precision",
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ]
            self.family_help.value = self.PSEUDO_HELP_WO_SOC

    def _update_family_link(self):
        library, accuracy = model.pseudos.library.split()
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

    def _update_pseudos(self):
        try:
            pseudo_family = self._get_pseudo_family()
            pseudos = pseudo_family.get_pseudos(structure=model.input_structure)
        except ValueError as exception:
            self._status_message.message = f"""
                <div class='alert alert-danger'>
                    ERROR: {exception!s}
                </div>
            """
            return

        model.pseudos.dictionary = {
            kind: pseudo.uuid for kind, pseudo in pseudos.items()
        }

    def _update_cutoffs(self):
        """Update wavefunction and density cutoffs from pseudo family."""
        if model.input_structure is None:
            model.pseudos.cutoffs = [[0.0], [0.0]]
            return

        try:
            pseudo_family = self._get_pseudo_family()
            current_unit = pseudo_family.get_cutoffs_unit()
            cutoff_dict = pseudo_family.get_cutoffs()
        except exceptions.NotExistent:
            self._status_message.message = f"""
                <div class='alert alert-danger'>
                    ERROR: required pseudo family `{model.pseudos.family}` is
                    not installed. Please use `aiida-pseudo install` to install
                    it."
                </div>
            """
        except ValueError as exception:
            self._status_message.message = f"""
                <div class='alert alert-danger'>
                    ERROR: failed to obtain recommended cutoffs for pseudos
                    `{pseudo_family}`: {exception}
                </div>
            """
        else:
            kind_names = model.input_structure.get_kind_names()

        ecutwfc_list = []
        ecutrho_list = []
        for kind in kind_names:
            cutoff = cutoff_dict.get(kind, {})
            ecutrho, ecutwfc = (
                U.Quantity(v, current_unit).to("Ry").to_tuple()[0]
                for v in cutoff.values()
            )
            ecutwfc_list.append(ecutwfc)
            ecutrho_list.append(ecutrho)

        model.pseudos.cutoffs = [ecutwfc_list, ecutrho_list]

    def _toggle_setter_widgets(self, change):
        if change["new"]:
            self.container.children += (self.setter_widget,)
        else:
            self.container.children = [*self.container.children][:1]

    def _build_setter_widgets(self):
        """Build the pseudo setter widgets."""

        self._unsubscribe_setter_links()

        children = []

        if model.input_structure is None:
            kinds = []
        else:
            self._update_pseudos()
            self._update_cutoffs()
            kinds = model.input_structure.kinds

        for index, kind in enumerate(kinds):
            symbol = kind.name
            upload_widget = PseudoUploadWidget(kind=symbol)
            pseudo_link = ipw.link(
                (model.pseudos, "dictionary"),
                (upload_widget, "pseudo"),
                [
                    lambda d, symbol=symbol: orm.load_node(d.get(symbol)),
                    lambda v, symbol=symbol: {
                        **model.pseudos.dictionary,
                        symbol: v.uuid,
                    },
                ],
            )
            cutoffs_link = ipw.dlink(
                (model.pseudos, "cutoffs"),
                (upload_widget, "cutoffs"),
                lambda c, i=index: [c[0][i], c[1][i]],
            )
            upload_widget.render()

            self.setter_widget_links.extend(
                [
                    pseudo_link,
                    cutoffs_link,
                    *upload_widget.links,
                ]
            )

            children.append(upload_widget)

        self.setter_widget.children = children

    def _get_pseudo_family(self):
        """Get the pseudo family from the database."""
        return (
            orm.QueryBuilder()
            .append(
                (
                    PseudoDojoFamily,
                    SsspFamily,
                    CutoffsPseudoPotentialFamily,
                ),
                filters={"label": model.pseudos.family},
            )
            .one()[0]
        )

    def _unsubscribe_setter_links(self):
        for link in self.setter_widget_links:
            link.unlink()
        self.setter_widget_links.clear()


class PseudoUploadWidget(ipw.HBox):
    """Class that allows to upload pseudopotential from user's computer."""

    pseudo = tl.Instance(UpfData, allow_none=True)
    cutoffs = tl.List(tl.Float, [])
    error_message = tl.Unicode(allow_none=True)

    def __init__(self, kind, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading pseudopotential uploader")],
            **kwargs,
        )

        self.kind = kind

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.file_upload = ipw.FileUpload(
            description="Upload",
            multiple=False,
        )
        self.pseudo_text = ipw.Text(description=self.kind)
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
            lambda p: p.filename if p else "",
        )

        cutoff_link = ipw.dlink(
            (self, "cutoffs"),
            (self.cutoff_message, "value"),
            lambda c: cutoffs_message_template.format(
                ecutwfc=c[0] if len(c) else "not set",
                ecutrho=c[1] if len(c) else "not set",
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
            element = "".join([i for i in self.kind if not i.isdigit()])
            if element != self.pseudo.element:
                self.error_message = f"""<div class="alert alert-danger"> ERROR: Element {self.kind} is not matched with the pseudo {self.pseudo.element}</div>"""
                self._reset()
            else:
                self.pseudo_text.value = filename

    def _reset(self):
        """Reset the widget to the initial state."""
        self.pseudo = None
        self.cutoffs = []
