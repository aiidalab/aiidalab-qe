import subprocess

import ipywidgets as ipw
import requests as req

from aiidalab_qe.app.utils import get_entry_items
from aiidalab_qe.common.widgets import RollingOutput
from aiidalab_widgets_base import LoadingWidget


class ArchiveImporter(ipw.VBox):
    GITHUB = "https://github.com"
    INFO_TEMPLATE = "{} <i class='fa fa-spinner fa-spin'></i>"
    DESCRIPTION_TEMPLATE = """
        <div class="alert alert-info" style="margin-bottom: 4px;">
            <h3>{header}</h3>
            <p>{content}</p>
        </div>
    """

    def __init__(
        self,
        repo: str,
        archive_list_url: str,
        archives_url: str,
        **kwargs,
    ):
        self.repo = repo
        self.archive_list_url = archive_list_url
        self.archives_url = archives_url
        self.archives: dict[str, dict[str, str]] = {}

        self.logger_placeholder = "Archive import output will be shown here."
        self.logger = RollingOutput()  # TODO use streaming output
        self.logger.value = self.logger_placeholder

        super().__init__(children=[LoadingWidget()], **kwargs)

    def render(self):
        self.selector = ipw.SelectMultiple(
            options=[],
            description="Examples:",
            rows=10,
            style={"description_width": "initial"},
            layout=ipw.Layout(width="auto"),
        )
        self.selector.observe(self._on_examples_selections, names="value")

        self.import_button = ipw.Button(
            description="Import",
            button_style="success",
            layout=ipw.Layout(width="fit-content"),
            icon="download",
        )
        self.import_button.on_click(self._on_import_click)

        self.info = ipw.HTML()

        accordion = None
        if self.logger:
            accordion = ipw.Accordion(children=[self.logger])
            accordion.set_title(0, "Archive import log")
            accordion.selected_index = None

        self.example_description = ipw.HTML()

        self.children = [
            ipw.HTML(f"""
                For questions regarding the examples, please open an issue in the
                <a href="{self.GITHUB}/{self.repo}" target="_blank">{self.repo}</a>
                repository.
            """),
            self.selector,
            ipw.HBox(
                children=[
                    self.import_button,
                    self.info,
                ],
                layout=ipw.Layout(margin="2px 0 4px 68px", grid_gap="4px"),
            ),
            self.example_description,
            accordion or ipw.Box(),
        ]

        self.selector.options = self._get_options()

    def _on_examples_selections(self, _) -> None:
        """Update the description of the selected example."""
        selected = self.selector.value
        if not selected:
            self.example_description.value = ""
            self.import_button.disabled = True
            return

        self.import_button.disabled = False

        if len(selected) > 1:
            description = self.DESCRIPTION_TEMPLATE.format(
                header="Multiple examples selected",
                content="To see a description, select only one example.",
            )
        else:
            archive = selected[0]
            metadata = self.archives[archive]
            description = self.DESCRIPTION_TEMPLATE.format(
                header=archive,
                content=metadata["description"] or metadata["label"],
            )

        self.example_description.value = description

    def _on_import_click(self, _):
        self.import_button.disabled = True
        self.logger.value = ""
        encountered_an_error = False
        for filename in self.selector.value:
            encountered_an_error = self._import_archive(filename)
        self.import_button.disabled = False
        if encountered_an_error:
            self.info.value = (
                "ERROR: some archives failed to import. See log for details"
            )

    def _import_archive(self, filename: str) -> bool:
        self.info.value = self.INFO_TEMPLATE.format(f"Importing {filename}")
        file_url = f"{self.archives_url}/{filename}"
        process = subprocess.Popen(
            ["verdi", "archive", "import", "-v", "critical", file_url],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        stdout, stderr = process.communicate()
        self._report(filename, stdout, stderr)
        return bool(stderr)

    def _report(self, filename: str, stdout: str, stderr: str):
        if stderr and "Success" not in stdout:
            self.logger.value += f"[ERROR] - importing {filename} failed\n\n{stderr}"
        else:
            self.info.value = ""
            self.logger.value += f"{stdout}"
        self.logger.value += f"\n{'#' * 80}\n\n"

    def _get_options(self) -> list[tuple[str, str]]:
        try:
            response: req.Response = req.get(self.archive_list_url)
            if not response.ok:
                self.info.value = "Failed to fetch archive list"
                return []
            self.archives = response.json()
            if not self.archives:
                self.info.value = "NOTE: Plugin does not yet provide examples"
            return [
                (
                    metadata["label"],
                    archive,
                )
                for archive, metadata in self.archives.items()
            ]
        except req.RequestException as e:
            self.info.value = f"Failed to fetch archive list: {e}"
            return []


class ExamplesImporter(ipw.Tab):
    def __init__(self, repo: str, tag: str):
        super().__init__()
        self.core_repo = repo
        self.core_tag = tag
        self.children = []
        self.titles = []
        self._load_tabs()

    def _archive_urls(self, repo: str, tag: str) -> tuple[str, str]:
        refs = f"refs/tags/{tag}"
        refs = "refs/heads/refactor"
        return (
            f"https://raw.githubusercontent.com/{repo}/{refs}/examples.json",
            f"https://github.com/{repo}/raw/{refs}/examples",
        )

    def _load_tabs(self):
        list_url, archives_url = self._archive_urls(
            repo=self.core_repo,
            tag=self.core_tag,
        )
        core_widget = ArchiveImporter(
            repo=self.core_repo,
            archive_list_url=list_url,
            archives_url=archives_url,
        )
        core_widget.render()
        self.children = [core_widget]
        self.set_title(0, "Core")

        entries: dict[str, dict] = get_entry_items(
            "aiidalab_qe.properties",
            "examples",
        )
        for i, (name, items) in enumerate(entries.items(), start=1):
            if not items:
                continue
            repo = items.get("repo", "")
            list_url, archives_url = self._archive_urls(
                repo=repo,
                tag=items.get("tag", ""),
            )

            plugin_widget = ArchiveImporter(
                repo=repo,
                archive_list_url=list_url,
                archives_url=archives_url,
            )
            plugin_widget.render()
            self.children += (plugin_widget,)
            self.set_title(i, items.get("title", name.capitalize().replace("_", " ")))
