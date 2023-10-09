# -*- coding: utf-8 -*-
from __future__ import annotations

import io
import re

import ipywidgets as ipw
import traitlets as tl
from aiida import orm
from aiida.common import exceptions
from aiida.plugins import DataFactory, GroupFactory
from aiida_quantumespresso.workflows.pw.base import PwBaseWorkChain
from aiidalab_widgets_base.utils import StatusHTML

from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS

UpfData = DataFactory("pseudo.upf")
SsspFamily = GroupFactory("pseudo.family.sssp")
PseudoDojoFamily = GroupFactory("pseudo.family.pseudo_dojo")
CutoffsPseudoPotentialFamily = GroupFactory("pseudo.family.cutoffs")


class PseudoFamilySelector(ipw.VBox):
    title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Accuracy and precision</h4></div>"""
    )

    description = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        The exchange-correlation functional and pseudopotential library is set by
        the <b>protocol</b> configured in the "Workflow" tab. Here you can
        override the defaults if desired.</div>""",
        layout=ipw.Layout(max_width="60%"),
    )

    pseudo_family_prompt = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; opacity:0.5;">
        <b><a href="https://www.materialscloud.org/discover/sssp/table/efficiency"
        target="_blank">Pseudopotential family</b></div>"""
    )
    pseudo_family_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px; opacity:0.5;">
        If you are unsure, select 'SSSP efficiency', which for
        most calculations will produce sufficiently accurate results at
        comparatively small computational costs. If your calculations require a
        higher accuracy, select 'SSSP accuracy' or 'PseudoDojo stringent', which will be computationally
        more expensive. SSSP is the standard solid-state pseudopotentials.
        The PseudoDojo used here has the SR relativistic type.</div>"""
    )

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
    disabled = tl.Bool()

    value = tl.Unicode(
        default_value=DEFAULT_PARAMETERS["advanced"]["pseudo_family"],
    )

    def __init__(self, **kwargs):
        # Enable manual setting of the pseudopotential family
        self._default_protocol = DEFAULT_PARAMETERS["workchain"]["protocol"]
        self.set_pseudo_family_prompt = ipw.HTML("<b>&nbsp;&nbsp;Override&nbsp;</b>")
        self.override_protocol_pseudo_family = ipw.Checkbox(
            description="",
            indent=False,
            value=False,
            layout=ipw.Layout(max_width="10%"),
        )
        self.set_pseudo_family_box = ipw.HBox(
            [self.set_pseudo_family_prompt, self.override_protocol_pseudo_family],
            layout=ipw.Layout(max_width="20%"),
        )
        self.show_ui = ipw.Valid(value=True)
        self.override_protocol_pseudo_family.observe(self.set_show_ui, "value")
        self.override_protocol_pseudo_family.observe(self.set_text_color, "value")
        self.override_protocol_pseudo_family.observe(self.set_value_trait, "value")

        # Choose the DFT functional
        self.dft_functional = ipw.Dropdown(
            options=["PBE", "PBEsol"],
            value=DEFAULT_PARAMETERS["advanced"]["pseudo_family"].split("/")[2],
            style={"description_width": "initial"},
        )
        self.dft_functional.observe(self.set_value_trait, "value")
        #
        pseudo_family_type = DEFAULT_PARAMETERS["advanced"]["pseudo_family"].split("/")[
            0
        ]
        if pseudo_family_type.upper() == "SSSP":
            pseudo_family_type += (
                " " + DEFAULT_PARAMETERS["advanced"]["pseudo_family"].split("/")[-1]
            )
        elif pseudo_family_type.upper() == "PSEUDODOJO":
            pseudo_family_type = "PseudoDojo " + pseudo_family_type.split("_")[-2]
        self.protocol_selection = ipw.ToggleButtons(
            options=[
                "SSSP efficiency",
                "SSSP precision",
                "PseudoDojo standard",
                "PseudoDojo stringent",
            ],
            value=pseudo_family_type,
            layout=ipw.Layout(max_width="80%"),
        )
        self.protocol_selection.observe(self.set_value_trait, "value")

        self.dft_functional_box = ipw.VBox(
            children=[
                self.dft_functional_prompt,
                self.dft_functional,
                self.dft_functional_help,
            ],
            layout=ipw.Layout(max_width="40%"),
        )
        self.pseudo_protocol_box = ipw.VBox(
            children=[
                self.pseudo_family_prompt,
                self.protocol_selection,
                self.pseudo_family_help,
            ],
            layout=ipw.Layout(max_width="60%"),
        )
        ipw.dlink((self.show_ui, "value"), (self.protocol_selection, "disabled"))
        ipw.dlink((self.show_ui, "value"), (self.dft_functional, "disabled"))

        self.set_value_trait()

        super().__init__(
            children=[
                self.title,
                ipw.HBox(
                    [self.description, self.set_pseudo_family_box],
                    layout=ipw.Layout(height="50px", justify_content="space-between"),
                ),
                ipw.HBox([self.dft_functional_box, self.pseudo_protocol_box]),
            ]
        )
        self.protocol = self._default_protocol

    def set_value_trait(self, _=None):
        family, protocol = self.protocol_selection.value.split()
        functional = self.dft_functional.value
        if family == "PseudoDojo":
            pseudo_family_string = f"PseudoDojo/0.4/{functional}/SR/{protocol}/upf"
        elif family == "SSSP":
            pseudo_family_string = f"SSSP/1.2/{functional}/{protocol}"
        else:
            raise ValueError(
                f"Unknown pseudo family {self.override_protocol_pseudo_family.value}"
            )
        self.value = pseudo_family_string

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
        self.protocol_selection.value = "SSSP efficiency"

    @tl.observe("protocol")
    def _protocol_changed(self, _):
        """Input protocol changed, update the widget values."""
        self._update_settings_from_protocol(self.protocol)

    def _update_settings_from_protocol(self, protocol):
        """Update the widget values from the given protocol, and trigger the callback."""
        pseudo_family = PwBaseWorkChain.get_protocol_inputs(protocol)["pseudo_family"]

        with self.hold_trait_notifications():
            # This changes will trigger callbacks
            self.set_from_pseudo_family(pseudo_family)

    def set_from_pseudo_family(self, pseudo_family):
        family, _, functional, accuracy = pseudo_family.split("/")
        protocol_selection = f"{family} {accuracy}"
        self.protocol_selection.value = protocol_selection
        self.dft_functional.value = functional


class PseudoSetter(ipw.VBox):
    structure = tl.Instance(klass=orm.StructureData, allow_none=True)
    pseudo_family = tl.Unicode(allow_none=True)

    # output pseudos
    pseudos = tl.Dict()

    # output cutoffs
    ecutwfc = tl.Float()
    ecutrho = tl.Float()

    _default_pseudo_setter_helper_text = """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; color: red;">
            Input structure is not set. Please set the structure first.
            </div>"""
    _update_pseudo_setter_helper_text = """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        The pseudopotential for each kind of atom in the structure can be set customly.
        The default pseudopotential and cutoffs are get from the pseudo family.
        The cutoffs used for the calculation are the maximum of the default from all pseudopotentials
        and can be set customly.
        </div>"""

    _cutoff_setter_helper_text = """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        Please set the cutoffs for the calculation. The default cutoffs are get from the pseudo family.
        </div>"""

    def __init__(
        self,
        structure: orm.StructureData | None = None,
        pseudo_family: str | None = None,
        **kwargs,
    ):
        self.pseudo_setting_widgets = ipw.VBox()
        self._status_message = StatusHTML(clear_after=20)

        # the initial cutoffs are set to 0
        self.ecutwfc = 0
        self.ecutrho = 0

        self.pseudo_setter_helper = ipw.HTML(self._default_pseudo_setter_helper_text)
        self.cutoff_setter_helper = ipw.HTML(self._cutoff_setter_helper_text)
        self.ecutwfc_setter = ipw.FloatText(
            description="Wavefunction cutoff (Ry)",
            style={"description_width": "initial"},
        )
        self.ecutrho_setter = ipw.FloatText(
            description="Charge density cutoff (Ry)",
            style={"description_width": "initial"},
        )
        self.ecutwfc_setter.observe(self._on_cutoff_change, names="value")
        self.ecutrho_setter.observe(self._on_cutoff_change, names="value")

        super().__init__(
            children=[
                self.pseudo_setter_helper,
                self.pseudo_setting_widgets,
                self.cutoff_setter_helper,
                ipw.HBox(
                    children=[
                        self.ecutwfc_setter,
                        self.ecutrho_setter,
                    ],
                ),
                self._status_message,
            ],
            **kwargs,
        )
        self._reset()
        with self.hold_trait_notifications():
            self.structure = structure
            self.pseudo_family = pseudo_family

    def _on_cutoff_change(self, _=None):
        self.ecutwfc = self.ecutwfc_setter.value
        self.ecutrho = self.ecutrho_setter.value

    def _reset_cutoff_widgets(self):
        self.ecutwfc_setter.value = 0
        self.ecutrho_setter.value = 0

    def _reset_traitlets(self):
        self.ecutwfc = 0
        self.ecutrho = 0
        self.pseudos = dict()

    def _reset(self):
        """Reset the pseudo setting widgets according to the structure
        by default the pseudos are get from the pseudo family
        """
        if self.structure is None:
            self._reset_cutoff_widgets()
            self._reset_traitlets()
            self.pseudo_setting_widgets.children = ()
            self.pseudo_setter_helper.value = self._default_pseudo_setter_helper_text
            return

        if self.pseudo_family is None:
            # this happened from the beginning when the widget is initialized
            # but also for the case when pseudo family is not provided which
            # won't happened for the real use but may happen for the test
            # so we still generate the pseudo setting widgets
            kinds = self.structure.get_kind_names()

            # Reset the traitlets, so the interface is clear setup
            self.pseudo_setting_widgets.children = ()
            self._reset_traitlets()

            # loop over the kinds and create the pseudo setting widget
            # (initialized with the pseudo from the family)
            for kind in kinds:
                pseudo_upload_widget = PseudoUploadWidget(kind=kind)

                # keep track of the changing of pseudo setting of each kind
                pseudo_upload_widget.observe(self._update_pseudos, ["pseudo"])
                self.pseudo_setting_widgets.children += (pseudo_upload_widget,)

            return

        try:
            pseudo_family = self._get_pseudos_family(self.pseudo_family)
        except exceptions.NotExistent as exception:
            self._status_message.message = (
                f"""<div class='alert alert-danger'> ERROR: {str(exception)}</div>"""
            )
            return

        try:
            pseudos = pseudo_family.get_pseudos(structure=self.structure)
            # get cutoffs dict of all elements
            cutoffs = self._get_cutoffs(pseudo_family)
        except ValueError as exception:
            self._status_message.message = f"""<div class='alert alert-danger'> ERROR: failed to obtain recommended cutoffs for pseudos `{pseudo_family}`: {exception}</div>"""
            return

        # success get family and cutoffs, set the traitlets accordingly
        # set the recommended cutoffs
        self.pseudos = {kind: pseudo.uuid for kind, pseudo in pseudos.items()}
        self.set_pseudos(self.pseudos, cutoffs)

    def _get_pseudos_family(self, pseudo_family: str) -> orm.Group:
        """Get the pseudo family from the database."""
        try:
            pseudo_set = (PseudoDojoFamily, SsspFamily, CutoffsPseudoPotentialFamily)
            pseudo_family = (
                orm.QueryBuilder()
                .append(pseudo_set, filters={"label": pseudo_family})
                .one()[0]
            )
        except exceptions.NotExistent as exception:
            raise exceptions.NotExistent(
                f"required pseudo family `{pseudo_family}` is not installed. Please use `aiida-pseudo install` to"
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

    def _create_pseudo_widget(self, kind):
        """The sigle line of pseudo setter widget"""
        return PseudoUploadWidget(kind=kind)

    @tl.observe("structure")
    def _structure_change(self, _):
        self._reset()
        self._update_pseudos()

    @tl.observe("pseudo_family")
    def _pseudo_family_change(self, _):
        self._reset()
        self._update_pseudos()

    def _update_pseudos(self, _=None):
        """Update the pseudos according to the pseudo setting widgets"""
        self._reset_cutoff_widgets()
        for w in self.pseudo_setting_widgets.children:
            if w.error_message is not None:
                self._status_message.message = w.error_message
                return

            if w.pseudo is not None:
                self.pseudos[w.kind] = w.pseudo.uuid
                self.pseudo_setter_helper.value = self._update_pseudo_setter_helper_text

                with self.hold_trait_notifications():
                    self.ecutwfc_setter.value = max(self.ecutwfc, w.ecutwfc)
                    self.ecutrho_setter.value = max(self.ecutrho, w.ecutrho)

    def set_pseudos(self, pseudos, cutoffs):
        # Reset the traitlets, so the interface is clear setup
        self.pseudo_setting_widgets.children = ()
        self._reset_traitlets()

        # loop over the kinds and create the pseudo setting widget
        # (initialized with the pseudo from the family)
        for kind in self.structure.kinds:
            element = kind.symbol
            pseudo = orm.load_node(pseudos.get(kind.name, None))
            _cutoffs = cutoffs.get(element, None)  # cutoffs for each element
            pseudo_upload_widget = PseudoUploadWidget(
                kind=kind.name, pseudo=pseudo, cutoffs=_cutoffs
            )

            # keep track of the changing of pseudo setting of each kind
            pseudo_upload_widget.observe(
                self._update_pseudos, ["pseudo", "ecutwfc", "ecutrho"]
            )
            self.pseudo_setting_widgets.children += (pseudo_upload_widget,)
        self._update_pseudos()


class PseudoUploadWidget(ipw.HBox):
    """Class that allows to upload pseudopotential from user's computer."""

    pseudo = tl.Instance(klass=UpfData, allow_none=True)
    kind = tl.Unicode()
    ecutwfc = tl.Float(allow_none=True)
    ecutrho = tl.Float(allow_none=True)
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
        self.kind = kind
        self.file_upload = ipw.FileUpload(
            description="Upload",
            multiple=False,
        )
        self.pseudo_text = ipw.Text(description=kind)
        self.file_upload.observe(self._on_file_upload, names="value")

        self._cutoff_message = ipw.HTML(
            self.cutoffs_message_template.format(ecutwfc=0, ecutrho=0)
        )

        if pseudo is not None:
            self.pseudo = pseudo
            self.pseudo_text.value = pseudo.filename

        self.error_message = None
        super().__init__(
            children=[
                self.pseudo_text,
                self.file_upload,
                self._cutoff_message,
            ],
            **kwargs,
        )
        # set the widget directly to trigger the traitlets set
        if cutoffs is not None:
            self.ecutwfc = cutoffs.get("cutoff_wfc", None)
            self.ecutrho = cutoffs.get("cutoff_rho", None)
            self._cutoff_message.value = self.cutoffs_message_template.format(
                ecutwfc=self.ecutwfc or "not set", ecutrho=self.ecutrho or "not set"
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
