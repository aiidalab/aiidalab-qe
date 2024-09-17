from __future__ import annotations

import io
import re

import ipywidgets as ipw
import traitlets as tl
from IPython.display import display

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


class PseudoFamilySelector(ipw.VBox):
    title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Accuracy and precision</h4></div>"""
    )
    PSEUDO_HELP_SOC = """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px; opacity:0.5;">
            Spin-orbit coupling (SOC) calculations are supported exclusively with PseudoDojo pseudopotentials.
            PseudoDojo offers these pseudopotentials in two versions: standard and stringent.
            Here, we utilize the FR (fully relativistic) type from PseudoDojo.
            Please ensure you choose appropriate cutoff values for your calculations.
        </div>"""

    PSEUDO_HELP_WO_SOC = """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px; opacity:0.5;">
        If you are unsure, select 'SSSP efficiency', which for
        most calculations will produce sufficiently accurate results at
        comparatively small computational costs. If your calculations require a
        higher accuracy, select 'SSSP accuracy' or 'PseudoDojo stringent', which will be computationally
        more expensive. SSSP is the standard solid-state pseudopotentials.
        The PseudoDojo used here has the SR relativistic type.</div>"""

    description = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        The exchange-correlation functional and pseudopotential library is set by
        the <b>protocol</b> configured in the "Workflow" tab. Here you can
        override the defaults if desired.</div>""",
        layout=ipw.Layout(max_width="60%"),
    )

    # XXX: the link is not correct after add pseudo dojo
    pseudo_family_prompt = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; opacity:0.5;">
        <b><a href="https://www.materialscloud.org/discover/sssp/table/efficiency"
        target="_blank">Pseudopotential family</b></div>"""
    )
    pseudo_family_help = ipw.HTML(PSEUDO_HELP_WO_SOC)

    dft_functional_prompt = ipw.HTML(
        """
        <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; opacity:0.5;"><b>
        Exchange-correlation  functional</b></div>"""
    )
    dft_functional_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px; opacity:0.5;">
        The exchange-correlation energy is calculated using this functional. We currently provide support for two
        well-established generalised gradient approximation (GGA) functionals:
        PBE and PBEsol.</div>"""
    )

    protocol = tl.Unicode(allow_none=True)
    spin_orbit = tl.Unicode()

    def __init__(self, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading pseudopotential selector")],
            **kwargs,
        )

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        # Enable manual setting of the pseudopotential family
        self.set_pseudo_family_prompt = ipw.HTML("<b>&nbsp;&nbsp;Override&nbsp;</b>")

        self.override = ipw.Checkbox(
            description="",
            indent=False,
            layout=ipw.Layout(max_width="10%"),
        )
        ipw.link(
            (model.pseudos, "override"),
            (self.override, "value"),
        )
        self.override.observe(self.set_show_ui, "value")
        self.override.observe(self.set_text_color, "value")
        self.override.observe(self.set_pseudo_family, "value")

        self.set_pseudo_family_box = ipw.HBox(
            [self.set_pseudo_family_prompt, self.override],
            layout=ipw.Layout(max_width="20%"),
        )

        self.show_ui = ipw.Valid(value=True)

        # the widget for DFT functional selection
        self.dft_functional = ipw.Dropdown(
            options=["PBE", "PBEsol"],
            style={"description_width": "initial"},
        )
        ipw.link(
            (model.pseudos, "dft_functional"),
            (self.dft_functional, "value"),
        )
        self.dft_functional.observe(self.set_pseudo_family, "value")

        self.library_selection = ipw.ToggleButtons(
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
            (self.library_selection, "value"),
        )
        self.library_selection.observe(self.set_pseudo_family, "value")

        self.dft_functional_box = ipw.VBox(
            children=[
                self.dft_functional_prompt,
                self.dft_functional,
                self.dft_functional_help,
            ],
            layout=ipw.Layout(max_width="40%"),
        )

        self.pseudo_setup_box = ipw.VBox(
            children=[
                self.pseudo_family_prompt,
                self.library_selection,
                self.pseudo_family_help,
            ],
            layout=ipw.Layout(max_width="60%"),
        )

        ipw.dlink((self.show_ui, "value"), (self.library_selection, "disabled"))
        ipw.dlink((self.show_ui, "value"), (self.dft_functional, "disabled"))

        self.children = [
            self.title,
            ipw.HBox(
                [self.description, self.set_pseudo_family_box],
                layout=ipw.Layout(height="50px", justify_content="space-between"),
            ),
            ipw.HBox([self.dft_functional_box, self.pseudo_setup_box]),
        ]

        ipw.dlink(
            (model, "protocol"),
            (self, "protocol"),
        )
        ipw.dlink(
            (model, "spin_orbit"),
            (self, "spin_orbit"),
        )

        self.rendered = True

    def set_pseudo_family(self, _=None):
        """The callback when the selection of pseudo family or dft functional is changed.
        Also triggered when the override checkbox is changed.
        This is the only method to set the value of the widget.
        """
        library, accuracy = model.pseudos.library.split()
        functional = model.pseudos.dft_functional
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
            # TODO previous version was referencing a nonexisting variable
            raise ValueError(f"Unknown pseudo family {library}")

        model.pseudos.family = pseudo_family_string

    def set_show_ui(self, change):
        self.show_ui.value = not change.new

    def set_text_color(self, change):
        opacity = 1.0 if change.new else 0.5

        for html in (
            self.pseudo_family_prompt,
            self.pseudo_family_help,
            self.dft_functional_help,
            self.dft_functional_prompt,
        ):
            old_opacity = re.match(
                r"[\s\S]+opacity:([\S]+);[\S\s]+", html.value
            ).groups()[0]
            html.value = html.value.replace(
                f"opacity:{old_opacity};", f"opacity:{opacity};"
            )

    def reset(self):
        """Reset the widget."""
        self._update_settings_from_protocol()

    @tl.observe("spin_orbit")
    def _update_library_selection(self, _):
        """Update the library selection according to the spin orbit value."""
        if model.spin_orbit == "soc":
            self.library_selection.options = [
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ]
            self.pseudo_family_help.value = self.PSEUDO_HELP_SOC
        else:
            self.library_selection.options = [
                "SSSP efficiency",
                "SSSP precision",
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ]
            self.pseudo_family_help.value = self.PSEUDO_HELP_WO_SOC

    @tl.observe("protocol")
    def _update_settings_from_protocol(self, _=None):
        """Update the widget values from the given protocol, and trigger the callback."""
        # FIXME: this rely on the aiida-quantumespresso, which is not ideal

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

        self.load_from_pseudo_family(pseudo_family)

    def load_from_pseudo_family(self, pseudo_family: PseudoFamily):
        """Reload the widget from the given pseudo family string."""
        with self.hold_trait_notifications():
            # will trigger the callback to set the value of widgets
            model.pseudos.library = f"{pseudo_family.library} {pseudo_family.accuracy}"
            model.pseudos.dft_functional = pseudo_family.functional


class PseudoSetter(ipw.VBox):
    input_structure = tl.Instance(klass=orm.StructureData, allow_none=True)
    pseudo_family = tl.Unicode(allow_none=True)
    override = tl.Bool(False)

    # output pseudos
    pseudos = tl.Dict()

    _update_pseudo_setter_helper_text = """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        The pseudopotential for each kind of atom in the structure can be set customly.
        The default pseudopotential and cutoffs are get from the pseudo family.
        The cutoffs used for the calculation are the maximum of the default from all pseudopotentials
        and can be set customly.
        </div>"""

    _cutoff_setter_helper_text = """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        Please set the cutoffs for the calculation. The default cutoffs are get from the pseudo family.
        </div>"""

    def __init__(self, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading pseudopotential setter")],
            **kwargs,
        )

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.setter_widget_out = ipw.Output()
        self.pseudo_setter_helper = ipw.HTML(self._update_pseudo_setter_helper_text)
        self.pseudo_setting_widgets = ipw.VBox()

        self._status_message = StatusHTML(clear_after=20)

        self.cutoff_setter_helper = ipw.HTML(self._cutoff_setter_helper_text)
        self.ecutwfc = ipw.FloatText(
            description="Wavefunction cutoff (Ry)",
            style={"description_width": "initial"},
        )
        ipw.link(
            (model.pseudos, "ecutwfc"),
            (self.ecutwfc, "value"),
        )
        self.ecutrho = ipw.FloatText(
            description="Charge density cutoff (Ry)",
            style={"description_width": "initial"},
        )
        ipw.link(
            (model.pseudos, "ecutrho"),
            (self.ecutrho, "value"),
        )

        self.children = [
            self.setter_widget_out,
            self.cutoff_setter_helper,
            ipw.HBox(
                children=[
                    self.ecutwfc,
                    self.ecutrho,
                ],
            ),
            self._status_message,
        ]

        # TODO why the hold? should we apply this to all the widgets?
        with self.hold_trait_notifications():
            ipw.dlink(
                (model, "input_structure"),
                (self, "input_structure"),
            )
            ipw.dlink(
                (model.pseudos, "family"),
                (self, "pseudo_family"),
            )
            ipw.dlink(
                (model.pseudos, "override"),
                (self, "override"),
            )

        self.rendered = True

    @tl.observe("input_structure", "pseudo_family", "override")
    def _render(self, _):
        self.setter_widget_out.clear_output()
        if model.input_structure and model.pseudos.override:
            with self.setter_widget_out:
                display(
                    ipw.VBox(
                        [
                            self.pseudo_setter_helper,
                            self.pseudo_setting_widgets,
                        ]
                    )
                )
        else:
            self.reset()

    def create_setter_widget(self, pseudos, cutoffs):
        if not model.input_structure:
            return None

        widgets_list = []
        for kind in model.input_structure.kinds:
            pseudo_upload_widget = PseudoUploadWidget(
                kind=kind.name,
                cutoffs=cutoffs.get(kind.symbol, None),
            )
            ipw.dlink(
                (model.pseudos, "dictionary"),
                (pseudo_upload_widget, "pseudo"),
                lambda d, kind=kind: d.get(
                    kind,
                    orm.load_node(pseudos.get(kind.name, None)),
                ),
            )
            pseudo_upload_widget.observe(
                lambda change, kind=kind: self._update_pseudos_dict(
                    {**model.pseudos.dictionary, kind: change["new"].uuid}
                ),
                "pseudo",
            )
            pseudo_upload_widget.render()
            widgets_list.append(pseudo_upload_widget)

        return ipw.VBox(widgets_list)

    def _update_pseudos_dict(self, pseudos):
        model.pseudos.dictionary = pseudos

    def _reset_cutoffs(self):
        """Reset cutoffs to 0"""
        model.pseudos.ecutwfc = 0.0
        model.pseudos.ecutrho = 0.0

    def reset(self):
        """Reset the pseudo setting widgets according to the structure
        by default the pseudos are get from the pseudo family
        """
        if self.input_structure is None:
            self._reset_cutoffs()
            model.pseudos.dictionary = {}
            return

        if model.pseudos.family is None:
            # this happened from the beginning when the widget is initialized
            # but also for the case when pseudo family is not provided which
            # won't happened for the real use but may happen for the test
            # so we still generate the pseudo setting widgets

            # Reset the traitlets, so the interface is clear setup
            model.pseudos.dictionary = {}

            # loop over the kinds and create the pseudo setting widget
            # (initialized with the pseudo from the family)
            widgets_list = []
            for kind in model.input_structure.get_kind_names():
                pseudo_upload_widget = PseudoUploadWidget(kind=kind)

                # keep track of the changing of pseudo setting of each kind
                pseudo_upload_widget.observe(self._update_widgets, "pseudo")
                widgets_list.append(pseudo_upload_widget)

            self.pseudo_setting_widgets.children = widgets_list

            return

        try:
            pseudo_family = self._get_pseudos_family()
        except exceptions.NotExistent as exception:
            self._status_message.message = (
                f"""<div class='alert alert-danger'> ERROR: {exception!s}</div>"""
            )
            return

        try:
            pseudos = pseudo_family.get_pseudos(structure=model.input_structure)
            cutoffs = self._get_cutoffs(pseudo_family)
        except ValueError as exception:
            self._status_message.message = f"""<div class='alert alert-danger'> ERROR: failed to obtain recommended cutoffs for pseudos `{pseudo_family}`: {exception}</div>"""
            return

        model.pseudos.dictionary = {
            kind: pseudo.uuid for kind, pseudo in pseudos.items()
        }
        self.set_pseudos(model.pseudos.dictionary, cutoffs)

    def _get_pseudos_family(self) -> orm.Group:
        """Get the pseudo family from the database."""
        try:
            pseudo_set = (PseudoDojoFamily, SsspFamily, CutoffsPseudoPotentialFamily)
            pseudo_family = (
                orm.QueryBuilder()
                .append(pseudo_set, filters={"label": model.pseudos.family})
                .one()[0]
            )
        except exceptions.NotExistent as exception:
            raise exceptions.NotExistent(
                f"required pseudo family `{model.pseudos.family}` is not installed. Please use `aiida-pseudo install` to"
                "install it."
            ) from exception

        return pseudo_family

    def _get_cutoffs(self, pseudo_family):
        """Get the cutoffs from the pseudo family."""
        from aiida_pseudo.common.units import U

        try:
            cutoffs = pseudo_family.get_cutoffs()
        except ValueError as exception:
            self._status_message.message = f"""<div class='alert alert-danger'> ERROR: failed to obtain recommended cutoffs for pseudos `{pseudo_family}`: {exception}</div>"""
            return

        current_unit = pseudo_family.get_cutoffs_unit()
        for element, cutoff in cutoffs.items():
            cutoffs[element] = {
                k: U.Quantity(v, current_unit).to("Ry").to_tuple()[0]
                for k, v in cutoff.items()
            }

        return cutoffs

    @tl.observe("input_structure", "pseudo_family")
    def _on_state_change(self, _):
        self.reset()
        self._update_widgets()

    @tl.observe("pseudos")
    def _update_widgets(self, _=None):
        """Update the pseudos according to the pseudo setting widgets"""
        self._reset_cutoffs()
        for w in self.pseudo_setting_widgets.children:
            if w.error_message is not None:
                self._status_message.message = w.error_message
                return

            if w.pseudo is not None:
                with self.hold_trait_notifications():
                    model.pseudos.ecutwfc = max(model.pseudos.ecutwfc, w.ecutwfc)
                    model.pseudos.ecutrho = max(model.pseudos.ecutrho, w.ecutrho)

    def set_pseudos(self, pseudos, cutoffs):
        self.pseudo_setting_widgets = self.create_setter_widget(pseudos, cutoffs)
        self._update_widgets()


class PseudoUploadWidget(ipw.HBox):
    """Class that allows to upload pseudopotential from user's computer."""

    pseudo = tl.Instance(UpfData, allow_none=True)
    error_message = tl.Unicode(allow_none=True)

    cutoffs_message_template = """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        The recommened ecutwfc: <b>{ecutwfc} Ry</b> &nbsp;
        for ecutrho: <b>{ecutrho} Ry</b>
        </div>"""

    def __init__(
        self,
        kind: str = "",
        pseudo: UpfData | None = None,
        cutoffs: dict | None = None,
        **kwargs,
    ):
        super().__init__(
            children=[LoadingWidget("Loading pseudopotential uploader")],
            **kwargs,
        )

        self.kind = kind
        self.cutoffs = cutoffs
        self.ecutwfc = None
        self.ecutrho = None
        if pseudo is not None:
            self.pseudo = pseudo

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

        self._cutoff_message = ipw.HTML(
            self.cutoffs_message_template.format(ecutwfc=0, ecutrho=0)
        )

        if self.pseudo is not None:
            self.pseudo_text.value = self.pseudo.filename

        self.error_message = None

        self.children = [
            self.pseudo_text,
            self.file_upload,
            self._cutoff_message,
        ]

        # set the widget directly to trigger the traitlets set
        if self.cutoffs is not None:
            self.ecutwfc = self.cutoffs.get("cutoff_wfc", None)
            self.ecutrho = self.cutoffs.get("cutoff_rho", None)
            self._cutoff_message.value = self.cutoffs_message_template.format(
                ecutwfc=self.ecutwfc or "not set",
                ecutrho=self.ecutrho or "not set",
            )

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
        self.ecutrho = None
        self.ecutwfc = None
