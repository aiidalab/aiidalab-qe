import base64
import pathlib
import tempfile
from threading import Thread

import ipywidgets as ipw


class DownloadDataWidget(ipw.VBox):
    def __init__(self, workchain_node):
        #
        self.download_archive_button = ipw.Button(
            description="Download AiiDA archive.aiida data",
            icon="download",
            button_style="primary",
            disabled=False,
            tooltip="Download the AiiDA archive of the simulation, ready to be shared or imported into another AiiDA profile",
            layout=ipw.Layout(width="100%"),
        )
        self.download_archive_button.on_click(self._download_data_thread)

        self.download_raw_button = ipw.Button(
            description="Download AiiDA raw data (zip format)",
            icon="download",
            button_style="primary",
            disabled=False,
            tooltip="Download the raw data of the simulation, organized in intuitive directory paths.",
            layout=ipw.Layout(width="100%"),
        )
        try:
            # check that we can import the ProcessDumper (not implemented in old AiiDA versions)
            # pre-commit: allow any unused imports in the next line
            from aiida.tools.dumping.processes import ProcessDumper  # noqa: F401

            self.download_raw_button.on_click(self._download_data_thread)
            self.dumper_is_available = True
        except Exception:
            self.dumper_is_available = False

        self.download_raw_button.disabled = not self.dumper_is_available

        self.node = workchain_node

        self._downloading_message = ipw.HTML()

        children = []

        if not self.dumper_is_available:
            children.append(
                ipw.HTML("""
                    <p style="color:red; line-height: 140%;">
                        The raw data download is not available because the AiiDA
                        version is too old.
                    </p>
                """),
            )

        children.extend(
            [
                ipw.HBox(children=[self.download_raw_button]),
                ipw.HBox(children=[self.download_archive_button]),
                self._downloading_message,
            ]
        )

        super().__init__(
            children=children,
        )

    def _download_data_thread(self, button_instance):
        thread = Thread(target=lambda: self._download_data(button_instance))
        thread.start()

    def _download_data(self, button_instance):
        """
        This method handles the download process when a download button is clicked.
        It updates the button's description to indicate that the download is in progress,
        determines whether to download the archive or raw data based on the button's description,
        generates the appropriate bitstream from the specified node, initiates the download
        with a filename based on the node's primary key, and then resets the button description
        to its original state.

        Args:
            button_instance (ipywidgets.Button): The button instance that was clicked.
        """
        if "archive" in button_instance.description:
            what = "archive"
            filename = f"export_qeapp_calculation_pk_{self.node.pk}.aiida"
        else:
            what = "raw"
            filename = f"export_{self.node.pk}_raw.zip"

        self._disable_buttons()
        self._show_downloading_message(what)
        data = self.produce_bitestream(self.node, what=what)
        self._download(payload=data, filename=filename)
        del data
        self._hide_downloading_message()
        self._enable_buttons()

    def _show_downloading_message(self, what):
        self._downloading_message.value = f"""
            <div style="display: flex; align-items: center; margin: 6px 0;">
                Creating {what} data to download
                <i class="fa fa-spinner fa-spin fa-2x fa-fw" style="margin-left: 4px;"></i>
            </div>
        """

    def _hide_downloading_message(self):
        self._downloading_message.value = ""

    def _disable_buttons(self):
        self.download_raw_button.disabled = True
        self.download_archive_button.disabled = True

    def _enable_buttons(self):
        if self.dumper_is_available:
            self.download_raw_button.disabled = False
        self.download_archive_button.disabled = False

    @staticmethod
    def _download(payload, filename):
        from IPython.display import Javascript, display

        javas = Javascript(
            f"""
            var link = document.createElement('a');
            link.href = 'data:application;base64,{payload}'
            link.download = '{filename}'
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """
        )
        display(javas)

    @staticmethod
    def produce_bitestream(node, what="archive"):
        """
        Produce a base64-encoded bitstream of the specified node data.

        Parameters:
        node (orm.Node): The AiiDA node to be processed.
        what (str): The type of data to produce. Options are "archive" or "raw".
                    Defaults to "archive".

        Returns:
        str: A base64-encoded string representing the requested data.

        Raises:
        KeyError: If the 'what' parameter is not "archive" or "raw".

        The function supports two modes:
        1. "archive": Creates an AiiDA archive of the node.
        2. "raw": Dumps the raw data of the process node into a zip file.

        NB: The function uses a temporary directory to store the data before converting it to a base64 string.
            Moreover, the node has to be reloaded because otherwise the SQLAlchemy will compleain on Db request
            not being in the same thread (the notebook session) as the original node.
        """
        from aiida import orm

        reloaded_node = orm.load_node(node.pk)
        with tempfile.TemporaryDirectory() as dirpath:
            if what == "archive":
                from aiida.tools.archive.create import create_archive

                path = pathlib.Path(dirpath) / "archive.aiida"
                create_archive(
                    entities=[reloaded_node],
                    filename=path,
                    call_calc_backward=False,
                    call_work_backward=False,
                    create_backward=False,
                )
                with open(path, "rb") as f:
                    zip_data = f.read()

                # Convert the ZIP data to base64 so it can be used as a payload in JavaScript
                bitestream = base64.b64encode(zip_data).decode()

            elif what == "raw":
                import shutil

                from aiida.tools.dumping.processes import ProcessDumper

                path = pathlib.Path(dirpath) / "raw_data"
                output_zip_path = pathlib.Path(dirpath) / "raw_data.zip"
                dumper = ProcessDumper()
                dumper.dump(process_node=reloaded_node, output_path=path)
                # writing files to a zipfile
                shutil.make_archive(pathlib.Path(dirpath) / "raw_data", "zip", path)

                with open(output_zip_path, "rb") as f:
                    raw_data = f.read()

                # Convert the raw_data to base64 so it can be used as a payload in JavaScript
                bitestream = base64.b64encode(raw_data).decode()

            else:
                raise KeyError("You should ask for `archive` or `raw` only!")

            return bitestream
