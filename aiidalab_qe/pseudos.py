# -*- coding: utf-8 -*-
import re
from subprocess import run
from threading import Thread

import ipywidgets as ipw
import traitlets
import yaml
from aiida import orm, plugins
from importlib_resources import files

from aiidalab_qe import parameters

SsspFamily = plugins.GroupFactory("pseudo.family.sssp")

DEFAULT_PARAMETERS = yaml.safe_load(files(parameters).joinpath("qeapp.yaml").open())


class Spinner(ipw.HTML):
    """Widget that shows a simple spinner if enabled."""

    enabled = traitlets.Bool()

    def __init__(self, spinner_style=None):
        self.spinner_style = f' style="{spinner_style}"' if spinner_style else ""
        super().__init__()

    @traitlets.default("enabled")
    def _default_enabled(self):  # pylint: disable=no-self-use
        return False

    @traitlets.observe("enabled")
    def _observe_enabled(self, change):
        """Show spinner if enabled, otherwise nothing."""
        if change["new"]:
            self.value = (
                f"""<i class="fa fa-spinner fa-pulse"{self.spinner_style}></i>"""
            )
        else:
            self.value = ""


class SSSPInstallWidget(ipw.HBox):

    installed = traitlets.Bool().tag(readonly=True)
    busy = traitlets.Bool().tag(readonly=True)

    base_url = "http://legacy-archive.materialscloud.org/file/2018.0001/v3"

    def __init__(self, **kwargs):

        self.install_button = ipw.Button(
            description="Install pseudos",
            button_style="warning",
            icon="cloud-download",
            tooltip="Download and install the SSSP pseudo potential families.",
            disabled=True,
            layout=ipw.Layout(width="140px"),
        )
        self.install_button.on_click(lambda _: self.download())

        self.spinner = Spinner()
        ipw.dlink((self, "busy"), (self.spinner, "enabled"))

        kwargs.setdefault("layout", ipw.Layout(width="180px"))

        super().__init__(children=[self.install_button, self.spinner], **kwargs)
        self._refresh_installed()

    def _refresh_installed(self):
        try:
            self.set_trait("busy", True)
            all_pseudo_families = (
                orm.QueryBuilder().append(SsspFamily, project="label").all(flat=True)
            )
            self.installed = (
                "SSSP/1.1/PBE/precision" in all_pseudo_families
                and "SSSP/1.1/PBE/efficiency" in all_pseudo_families
                and "SSSP/1.1/PBEsol/precision" in all_pseudo_families
                and "SSSP/1.1/PBEsol/efficiency" in all_pseudo_families
            )
        finally:
            self.set_trait("busy", False)

    @traitlets.default("busy")
    def _default_busy(self):
        return False

    @traitlets.observe("busy")
    def _observe_busy(self, change):
        self.install_button.disabled = change["new"]

    @traitlets.observe("installed")
    def _observe_installed(self, change):
        self.install_button.layout.visibility = "hidden" if change["new"] else "visible"

    def _download(self):
        try:
            self.set_trait("busy", True)

            functionals = ["PBE", "PBEsol"]
            protocols = ["efficiency", "precision"]

            for func in functionals:
                for prot in protocols:
                    run(["aiida-pseudo", "install", "sssp", "-x", func, "-p", prot])
            self._refresh_installed()
        finally:
            self.set_trait("busy", False)

    def download(self):
        thread = Thread(target=self._download)
        thread.start()


class PseudoFamilySelector(ipw.VBox):

    title = ipw.HTML(
        """<div style="padding-top: 0px; padding-bottom: 10px">
        <h4>Accuracy and precision</h4></div>"""
    )

    description = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px">
        The exchange-correlation function and pseudopotential library is set by
        the <b>protocol</b> configured in the "WorkChain" tab. Here you can
        override the defaults if desired.</div>""",
        layout=ipw.Layout(max_width="60%"),
    )

    pseudo_family_prompt = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; opacity:0.5;">
        <b><a href="https://www.materialscloud.org/discover/sssp/table/precision"
        target="_blank">Standard Solid-state pseudopotential</a> (SSSP) protocol</b></div>"""
    )
    pseudo_family_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 0px; opacity:0.5;">
        If you are unsure what to choose, select 'SSSP efficiency', which for
        most calculations will produce sufficiently accurate results at
        comparatively small computational cost. If your calculation requires a
        higher accuracy, select 'SSSP accuracy', which will be computationally
        more expensive, but will produce even more accurate results.</div>"""
    )

    dft_functional_prompt = ipw.HTML(
        """
        <div style="line-height: 140%; padding-top: 0px; padding-bottom: 10px; opacity:0.5;"><b>
        Exchange-correlation  functional</b></div>"""
    )
    dft_functional_help = ipw.HTML(
        """<div style="line-height: 140%; padding-top: 10px; padding-bottom: 10px; opacity:0.5;">
        The exchange-correlation energy is calculated based on the charge
        density using this functional. We currently provide support for two
        well-established generalised gradient approximation (GGA) functionals:
        PBE and PBEsol.</div>"""
    )
    installed = traitlets.Bool()
    disable = traitlets.Bool()

    value = traitlets.Unicode(
        default_value=DEFAULT_PARAMETERS["pseudo_family"],
    )
    FUNCTIONAL_DEFAULT = DEFAULT_PARAMETERS["pseudo_family"].split("/")[2]
    PROTOCOL_DEFAULT = DEFAULT_PARAMETERS["pseudo_family"].split("/")[3]

    def __init__(self, **kwargs):

        # Enable manual setting of the pseudopotential family
        self.set_pseudo_family_prompt = ipw.HTML("<b>&nbsp;&nbsp;Override&nbsp;</b>")
        self.set_pseudo_family = ipw.Checkbox(
            description="",
            indent=False,
            value=False,
            layout=ipw.Layout(max_width="10%"),
        )
        self.set_pseudo_family_box = ipw.HBox(
            [self.set_pseudo_family_prompt, self.set_pseudo_family],
            layout=ipw.Layout(max_width="20%"),
        )
        self.show_ui = ipw.Valid(value=True)
        self.set_pseudo_family.observe(self.set_show_ui, "value")
        self.set_pseudo_family.observe(self.set_text_color, "value")

        # Choose the DFT functional
        self.dft_functional = ipw.Dropdown(
            options=["PBE", "PBEsol"],
            value=self.FUNCTIONAL_DEFAULT,
            style={"description_width": "initial"},
        )
        self.dft_functional.observe(self.set_value_trait, "value")

        self.protocol_selection = ipw.ToggleButtons(
            options=["efficiency", "precision"],
            value=self.PROTOCOL_DEFAULT,
            layout=ipw.Layout(max_width="80%"),
        )
        self.protocol_selection.observe(self.set_value_trait, "value")

        # Setup pseudofamily potential selection group:
        self.sssp_install_widget = SSSPInstallWidget()
        ipw.dlink((self.sssp_install_widget, "installed"), (self, "installed"))

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

    def set_value_trait(self, _=None):
        self.value = (
            f"SSSP/1.1/{self.dft_functional.value}/{self.protocol_selection.value}"
        )

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
        self.protocol_selection.value = "efficiency"
