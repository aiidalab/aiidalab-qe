from __future__ import annotations

import io

import ipywidgets as ipw

from aiidalab_qe.utils import generate_alert
from aiidalab_widgets_base import LoadingWidget
from aiidalab_widgets_base.utils import StatusHTML

from ..utils import (
    UpfData,
    get_pseudo_by_filename,
    get_pseudo_by_md5,
)
from .model import PseudoPotentialUploaderModel


class PseudoPotentialUploader(ipw.VBox):
    """Class that allows to upload pseudopotential from user's computer."""

    def __init__(self, model: PseudoPotentialUploaderModel, **kwargs):
        super().__init__(
            children=[LoadingWidget("Loading pseudopotential uploader")],
            **kwargs,
        )
        self._model = model

        self._model.observe(
            self._on_pseudo_change,
            "pseudo",
        )

        self.rendered = False

    def render(self):
        if self.rendered:
            return

        self.pseudo_filename = ipw.Text(
            description=self._model.kind_name,
            style={"description_width": "50px"},
            disabled=True,
        )
        filename_link = ipw.dlink(
            (self._model, "pseudo"),
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
        self.file_upload.observe(
            self._on_file_upload,
            "value",
        )

        self.pseudo_info = ipw.HTML()
        info_link = ipw.dlink(
            (self._model, "info"),
            (self.pseudo_info, "value"),
        )

        self.message_box = StatusHTML(clear_after=5)
        message_link = ipw.link(
            (self._model, "message"),
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
        self._model.update_pseudo_info()

    def _on_pseudo_change(self, _):
        self._model.populate_extras()

    def _on_file_upload(self, change=None):
        if not change or not change["new"]:
            return

        filename, item = next(iter(change["new"].items()))
        content = item["content"]

        try:
            uploaded_pseudo = UpfData(io.BytesIO(content), filename=filename)
        except Exception:
            self._model.message = generate_alert(
                alert_type="danger",
                message=f"{filename} is not a valid UPF file",
            )
            self._reset_uploader()
            return

        # Wrong element
        if uploaded_pseudo.element != self._model.kind_symbol:
            self._model.message = generate_alert(
                alert_type="danger",
                message=f"Pseudo element {uploaded_pseudo.element} does not match {self._model.kind_symbol}",
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
            self._model.message = generate_alert(
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

        self._model.pseudo = uploaded_pseudo
        self._model.message = message
        self._model.uploaded = True
        self._reset_uploader()

    def _reset_uploader(self):
        self.file_upload.metadata = []
        self.file_upload.data = []
        self.file_upload._counter = 0
