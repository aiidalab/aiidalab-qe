# -*- coding: utf-8 -*-
import ipywidgets as ipw
from threading import Thread
from subprocess import run
from aiida import orm, plugins

import traitlets

SsspFamily = plugins.GroupFactory("pseudo.family.sssp")


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

    pseudo_family_prompt = ipw.HTML(
        'Select the <a href="https://www.materialscloud.org/discover/sssp/table/precision" '
        'target="_blank">pseudopotential library</a> for the calculation.'
    )

    pseudo_family_help = ipw.HTML(
        """
        <div style="line-height:120%;">If you are unsure what to choose, select 'SSSP efficiency', which for most
        calculations will produce sufficiently accurate results at comparatively small computational cost. If
        your calculation requires a higher accuracy, select 'SSSP accuracy', which will be computationally more
        expensive, but will produce even more accurate results.</div>"""
    )

    installed = traitlets.Bool()
    disabled = traitlets.Bool()

    value = traitlets.Unicode(allow_none=True)

    def __init__(self, **kwargs):
        self.protocol_selection = ipw.ToggleButtons(options=["efficiency", "precision"])
        self.protocol_selection.observe(self.set_value_trait, "value")

        # Setup pseudofamily potential selection group:
        self.sssp_install_widget = SSSPInstallWidget()
        ipw.dlink((self.sssp_install_widget, "installed"), (self, "installed"))

        # DFT functional.
        self.dft_functional = ipw.Dropdown(
            options=["PBE", "PBEsol"],
            value="PBE",
            style={"description_width": "initial"},
        )
        self.dft_functional.observe(self.set_value_trait, "value")

        self.dft_functional_box = ipw.VBox(
            children=[
                ipw.HBox(
                    children=[
                        ipw.HTML(
                            "<b>DFT functional</b>", layout=ipw.Layout(flex="1 1 auto")
                        ),
                        self.dft_functional,
                    ]
                ),
                ipw.HTML("""Some explanation for the functional.... """),
            ],
            layout=ipw.Layout(max_width="600px"),
        )

        # Pseudo potential family selection.
        self.chose_automatically = ipw.Checkbox(
            description="Chose pseudo automatically",
            indent=False,
            value=True,
        )
        ipw.dlink(
            (self.chose_automatically, "value"), (self.protocol_selection, "disabled")
        )
        ipw.dlink(
            (self.chose_automatically, "value"), (self.dft_functional, "disabled")
        )
        self.chose_automatically.observe(self.set_value_trait, "value")
        self.set_value_trait()

        super().__init__(
            children=[
                self.chose_automatically,
                self.dft_functional_box,
                self.pseudo_family_prompt,
                ipw.HBox([self.protocol_selection, self.sssp_install_widget]),
                self.pseudo_family_help,
            ]
        )

    def set_value_trait(self, _=None):
        self.value = (
            None
            if self.chose_automatically.value
            else f"SSSP/1.1/{self.dft_functional.value}/{self.protocol_selection.value}"
        )

    def reset(self):
        self.protocol_selection.value = "efficiency"
