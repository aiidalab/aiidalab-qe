import ipywidgets as ipw


def get_start_widget(appbase, jupbase, notebase):  # noqa: ARG001
    return ipw.HTML(f"""
        <div class="app-container">
            <a
                class="logo"
                href="{appbase}/qe.ipynb"
                target="_blank"
            >
                <img src="{appbase}/src/aiidalab_qe/app/static/images/logo.png" />
            </a>
            <div class="features">
                <a
                    class="feature"
                    href="{appbase}/qe.ipynb"
                    target="_blank">
                    <img
                        class="feature-logo"
                        src="{appbase}/miscellaneous/logos/workbench.png"
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
                        src="{appbase}/miscellaneous/logos/history.png"
                        alt="Calculation history"
                    />
                    <div class="feature-label">Calculation History</div>
                </a>
                <a
                    class="feature"
                    href="{appbase}/plugin_list.ipynb"
                    target="_blank">
                    <img
                        class="feature-logo"
                        src="{appbase}/miscellaneous/logos/plugins.png"
                        alt="Plugin store"
                    />
                    <div class="feature-label">Plugin store</div>
                </a>
                <a
                    class="feature"
                    href="https://www.quantum-espresso.org/"
                    target="_blank">
                    <img
                        class="feature-logo"
                        src="{appbase}/miscellaneous/logos/qe-logo.png"
                        alt="Quantum ESPRESSO"
                    />
                    <div class="feature-label">Quantum ESPRESSO</div>
                </a>
            </div>
        </div>
    """)
