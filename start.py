import base64
import mimetypes
from pathlib import Path

import ipywidgets as ipw


_APP_DIR = Path(__file__).parent


def _image_data_uri(relative_path):
    path = _APP_DIR / relative_path
    mime_type = mimetypes.guess_type(path.name)[0] or "image/png"
    encoded = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime_type};base64,{encoded}"


def get_start_widget(appbase, jupbase, notebase):  # noqa: ARG001
    logo = _image_data_uri("src/aiidalab_qe/app/static/images/logo.png")
    workbench = _image_data_uri("miscellaneous/logos/workbench.png")
    history = _image_data_uri("miscellaneous/logos/history.png")
    plugins = _image_data_uri("miscellaneous/logos/plugins.png")
    download = _image_data_uri("miscellaneous/logos/download.png")
    qe_logo = _image_data_uri("miscellaneous/logos/qe-logo.png")

    return ipw.HTML(f"""
        <style>
            details {{
                border: 1px solid #aaaaaa;
                border-radius: 4px;
                padding: 1rem 2rem 0.25rem;
                margin: 1rem 2px 2rem;
            }}
            summary {{
                font-weight: bold;
                margin: -0.5em -0.5em 0;
                padding: 0.5em;
                cursor: pointer;
            }}
            summary::before {{
                content: "▶";
            }}
            details[open] summary::before {{
                content: "▼";
            }}
        </style>
        <div class="app-container">
            <a
                class="logo"
                href="{appbase}/qe.ipynb"
                target="_blank"
            >
                <img src="{logo}" />
            </a>
            <div class="features">
                <a
                    class="feature"
                    href="{appbase}/qe.ipynb"
                    target="_blank">
                    <img
                        class="feature-logo"
                        src="{workbench}"
                        alt="New calculation"
                    />
                    <div class="feature-label">New calculation</div>
                </a>
                <a
                    class="feature"
                    href="{appbase}/calculation_history.ipynb"
                    target="_blank">
                    <img
                        class="feature-logo"
                        src="{history}"
                        alt="Calculation history"
                    />
                    <div class="feature-label">Calculation history</div>
                </a>
                <a
                    class="feature"
                    href="{appbase}/plugin_manager.ipynb"
                    target="_blank">
                    <img
                        class="feature-logo"
                        src="{plugins}"
                        alt="Plugin store"
                    />
                    <div class="feature-label">Plugin store</div>
                </a>
                <a
                    class="feature"
                    href="{appbase}/examples.ipynb"
                    target="_blank">
                    <img
                        class="feature-logo"
                        src="{download}"
                        alt="Download examples"
                    />
                    <div class="feature-label">Download examples</div>
                </a>
                <a
                    class="feature"
                    href="https://www.quantum-espresso.org/"
                    target="_blank">
                    <img
                        class="feature-logo"
                        src="{qe_logo}"
                        alt="Quantum ESPRESSO"
                    />
                    <div class="feature-label">Quantum ESPRESSO</div>
                </a>
            </div>

        </div>
    """)
