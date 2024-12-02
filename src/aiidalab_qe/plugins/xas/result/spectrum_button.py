import base64
import hashlib
from typing import Callable

import ipywidgets as ipw
from IPython.display import HTML, display


class SpectrumDownloadButton(ipw.Button):
    """Download button with dynamic content
    The content is generated using a callback when the button is clicked.
    Modified from responses to https://stackoverflow.com/questions/61708701/how-to-download-a-file-using-ipywidget-button#62641240
    """

    def __init__(self, filename: str, contents: Callable[[], str], **kwargs):
        super().__init__(**kwargs)
        self.filename = filename
        self.contents = contents
        self.on_click(self.__on_click)

    def __on_click(self, _):
        if self.contents is None:
            return

        contents: bytes = self.contents().encode("utf-8")
        b64 = base64.b64encode(contents)
        payload = b64.decode()
        digest = hashlib.md5(contents).hexdigest()  # bypass browser cache
        link_id = f"dl_{digest}"

        display(
            HTML(
                f"""
            <html>
            <body>
            <a id="{link_id}" download="{self.filename}" href="data:text/csv;base64,{payload}" download>
            </a>
            <script>
            (function download() {{
            document.getElementById('{id}').click();
            }})()
            </script>
            </body>
            </html>
        """
            )
        )
