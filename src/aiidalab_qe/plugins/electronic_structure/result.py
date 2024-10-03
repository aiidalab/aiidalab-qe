"""Electronic structure results view widgets"""

import ipywidgets as ipw
from IPython.display import Javascript, display

from aiidalab_qe.common.bandpdoswidget import BandPdosWidget
from aiidalab_qe.common.panel import ResultPanel


class Result(ResultPanel):
    title = "Electronic Structure"
    workchain_labels = ["bands", "pdos"]

    def __init__(self, node=None, **kwargs):
        super().__init__(node=node, **kwargs)

        # Create an HTML element for the warning message
        self.message_element = ipw.HTML(
            """<div id="window-size-message" style="color: red; display: none;"></div>"""
        )
        self.children = [self.message_element]

    def _update_view(self):
        """Update the view of the widget."""
        #
        try:
            pdos_node = self.node.outputs.pdos
        except AttributeError:
            pdos_node = None

        try:
            bands_node = self.node.outputs.bands
        except AttributeError:
            bands_node = None
        _bands_dos_widget = BandPdosWidget(bands=bands_node, pdos=pdos_node)

        # Create a unique ID for this instance
        unique_id = id(_bands_dos_widget)  # Use the object's id as a unique identifier

        # Create the message HTML element
        self.message_element = ipw.HTML(
            f"<div id='window-size-message-{unique_id}' style='color: red; display: none;'></div>"
        )

        # Display the message element first
        display(self.message_element)

        # Inject JavaScript to check window size and show the message
        js_code = f"""
        (function() {{
            function checkWindowSize() {{
                var width = window.innerWidth;
                var messageElement = document.getElementById('window-size-message-{unique_id}');

                // Check if messageElement exists
                if (!messageElement) {{
                    console.error("Message element not found.");
                    return;
                }}

                // Assuming only one plotly graph on the page; adjust if necessary
                var plotElement = document.getElementsByClassName('plotly-graph-div')[0];

                if (width < 1080) {{  // Adjust this value as needed
                    messageElement.innerHTML = "NOTE: The Electronic Plot window is too narrow - enlarge it to see also the legend and toggle the plot lines.";
                    messageElement.style.display = "block";

                }} else {{
                    messageElement.style.display = "none"; // Hide message if the window size is sufficient
                }}
            }}

            // Check the window size on load
            window.onload = checkWindowSize;

            // Add an event listener to check the size on window resize
            window.onresize = checkWindowSize;

            // Initial check to display message if necessary
            checkWindowSize();
        }})();
        """

        # Display the JavaScript after the message element
        display(Javascript(js_code))

        # update the electronic structure tab
        self.children = [_bands_dos_widget, self.message_element]
