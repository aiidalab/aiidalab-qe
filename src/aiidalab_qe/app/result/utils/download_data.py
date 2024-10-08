import base64
import pathlib
import tempfile

import ipywidgets as ipw


class DownloadDataWidget(ipw.VBox):
    def __init__(self, qeapp_node):
        self.download_archive_button = ipw.Button(
            description="Download AiiDA archive zip data",
            icon="download",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )
        self.download_archive_button.on_click(self.download_archive_data)

        self.download_raw_button = ipw.Button(
            description="Download AiiDA raw zip data",
            icon="download",
            button_style="primary",
            disabled=False,
            layout=ipw.Layout(width="auto"),
        )
        self.download_raw_button.on_click(self.download_raw_data)

        self.node = qeapp_node

        super().__init__(
            children=[
                self.download_archive_button,
                self.download_raw_button,
            ],
        )

    def download_archive_data(self, _=None):
        """
        Download both the phonopy.yaml and fc.hdf5 files.
        """
        self.download_archive_button.description += "... Downloading now..."
        archive_data = self.produce_bitestream(self.node, what="archive")
        self._download(payload=archive_data, filename=f"export_{self.node.pk}.aiida")
        del archive_data
        self.download_archive_button.description = (
            self.download_archive_button.description.replace(
                "... Downloading now...", ""
            )
        )

    def download_raw_data(self, _=None):
        """
        Download both the phonopy.yaml and fc.hdf5 files.
        """
        self.download_raw_button.description += "... Downloading now..."
        raw_data = self.produce_bitestream(self.node, what="raw")
        self._download(payload=raw_data, filename=f"export_{self.node.pk}_raw.zip")
        del raw_data
        self.download_raw_button.description = (
            self.download_raw_button.description.replace("... Downloading now...", "")
        )

    @staticmethod
    def _download(payload, filename):
        from IPython.display import Javascript, display

        javas = Javascript(
            f"""
            var link = document.createElement('a');
            link.href = 'data:application/octet-stream;base64,{payload}'
            link.download = "{filename}"
            document.body.appendChild(link);
            link.click();
            document.body.removeChild(link);
            """
        )
        display(javas)

    @staticmethod
    def produce_bitestream(node, what="archive"):
        with tempfile.TemporaryDirectory() as dirpath:
            if what == "archive":
                from aiida.tools.archive.create import create_archive

                path = pathlib.Path(dirpath) / "archive.aiida"
                create_archive(entities=[node], filename=path)
                with open(path, "rb") as f:
                    zip_data = f.read()

                # Convert the ZIP data to base64 so it can be used as a payload in JavaScript
                bitestream = base64.b64encode(zip_data).decode("utf-8")

            elif what == "raw":
                import shutil

                from aiida.tools.dumping.processes import ProcessDumper

                path = pathlib.Path(dirpath) / "raw_data"
                output_zip_path = pathlib.Path(dirpath) / "raw_data.zip"
                ProcessDumper().dump(process_node=node, output_path=path)
                # writing files to a zipfile
                shutil.make_archive(pathlib.Path(dirpath) / "raw_data", "zip", path)

                with open(output_zip_path, "rb") as f:
                    raw_data = f.read()

                # Convert the raw_data to base64 so it can be used as a payload in JavaScript
                bitestream = base64.b64encode(raw_data).decode("utf-8")

            else:
                raise KeyError("You should ask for `archive` or `raw` only!")

            return bitestream
