from datetime import datetime, timezone

import ipywidgets as ipw
from table_widget import TableWidget

from aiida.orm import QueryBuilder, load_node
from aiidalab_qe.common.process import STATE_ICONS
from aiidalab_qe.common.widgets import LoadingWidget


def determine_state_icon(row):
    """Attach an icon to the displayed job state."""
    state = row["state"].lower()
    if state == "finished" and row.get("exit_status", 0) != 0:
        return f"Finished{STATE_ICONS['excepted']}"
    return f"{state.capitalize()}{STATE_ICONS.get(state, '')}"


COLUMNS = {
    "id": {"headerName": "ID 🔗", "dataType": "link", "editable": False},
    "creation_time_absolute": {
        "headerName": "Creation time ⏰ (absolute)",
        "type": "date",
        "width": 100,
        "editable": False,
    },
    "creation_time_relative": {
        "headerName": "Creation time ⏰ (relative)",
        "width": 100,
        "editable": False,
    },
    "structure": {"headerName": "Structure", "editable": False},
    "state": {"headerName": "State 🟢", "editable": False},
    "status": {"headerName": "Status", "editable": False, "hide": True},
    "exit_status": {
        "headerName": "Exit status",
        "type": "number",
        "editable": False,
        "hide": True,
    },
    "exit_message": {"headerName": "Exit message", "editable": False},
    "label": {"headerName": "Label", "width": 300, "editable": True},
    "description": {
        "headerName": "Description",
        "width": 300,
        "editable": True,
        "hide": True,
    },
    "relax_type": {"headerName": "Relax type", "editable": False, "hide": True},
    "delete": {"headerName": "Delete", "dataType": "link", "editable": False},
    "download": {"headerName": "Download", "dataType": "link", "editable": False},
    "uuid": {"headerName": "UUID", "editable": False, "hide": True},
    "properties": {"headerName": "Properties", "editable": False, "hide": True},
}


class CalculationHistory:
    def __init__(self):
        self.main = ipw.VBox(children=[LoadingWidget("Loading the table...")])

        self.table = TableWidget(style={"margin-top": "20px"})

        def on_row_update(change):
            # When the user updates 'label' or 'description' in the table,
            # reflect these changes in the corresponding AiiDA node.
            node = load_node(change["new"]["uuid"])
            node.label = change["new"]["label"]
            node.description = change["new"]["description"]

        self.table.observe(on_row_update, "updatedRow")

        # This will hold the raw data (list of dicts) from the database
        self.data = []
        # This will hold the currently displayed data (filtered, with toggles applied, etc.)
        self.display_data = []

    def load_table(self):
        """Populate the table after initialization."""
        self.data = self.load_data()
        self.setup_widgets()

    def load_data(self):
        """Fetch the QeAppWorkChain results using the QueryBuilder and
        return a list of dictionaries with all required columns.

        """
        from aiidalab_qe.workflows import QeAppWorkChain

        projections = [
            "id",
            "uuid",
            "extras.structure",
            "ctime",
            "attributes.process_state",
            "attributes.process_status",
            "attributes.exit_status",
            "attributes.exit_message",
            "label",
            "description",
            "extras.workchain.relax_type",
            "extras.workchain.properties",
        ]

        qb = QueryBuilder()
        qb.append(QeAppWorkChain, project=projections, tag="process")
        qb.order_by({"process": {"ctime": "desc"}})
        results = qb.all()

        data = []
        if not results:
            return data

        now = datetime.now(timezone.utc)

        for row in results:
            (
                pk,
                uuid,
                structure,
                creation_time,
                state,
                status,
                exit_status,
                exit_message,
                label,
                description,
                relax_type,
                properties,
            ) = row

            creation_time_str = (
                creation_time.strftime("%Y-%m-%d %H:%M:%S") if creation_time else ""
            )
            if creation_time:
                days_ago = (now - creation_time).days
                creation_time_rel = f"{days_ago}D ago"
            else:
                creation_time_rel = "N/A"

            # Transform "waiting" to "running" for readbility
            if state == "waiting":
                state = "running"

            # Prepare link-based values
            pk_with_link = f'<a href="./qe.ipynb?pk={pk}" target="_blank">{pk}</a>'
            uuid_with_link = (
                f'<a href="./qe.ipynb?pk={pk}" target="_blank">{uuid[:8]}</a>'
            )
            delete_link = f'<a href="./delete.ipynb?pk={pk}" target="_blank">Delete</a>'
            download_link = (
                f'<a href="./download.ipynb?pk={pk}" target="_blank">Download</a>'
            )

            # Make sure properties is a list (avoid None)
            properties = properties if properties is not None else []

            data.append(
                {
                    "pk_with_link": pk_with_link,
                    "uuid_with_link": uuid_with_link,
                    "creation_time_absolute": creation_time_str,
                    "creation_time_relative": creation_time_rel,
                    "structure": structure,
                    "state": state,
                    "status": status,
                    "exit_status": exit_status,
                    "exit_message": exit_message,
                    "label": label,
                    "description": description,
                    "relax_type": relax_type,
                    "uuid": uuid,
                    "delete": delete_link,
                    "download": download_link,
                    "properties": properties,
                    "creation_time": creation_time,
                }
            )

        return data

    def setup_widgets(self):
        """Create widgets for filtering, toggles for display, etc."""
        # Gather unique properties
        all_properties = set()
        for row in self.data:
            for prop in row["properties"]:
                if prop is not None:
                    all_properties.add(prop)

        # Build a set of checkboxes for properties
        property_checkboxes = [
            ipw.Checkbox(
                value=False,
                description=prop,
                indent=False,
                layout=ipw.Layout(description_width="initial"),
            )
            for prop in sorted(all_properties)
        ]

        self.properties_box = ipw.HBox(property_checkboxes)
        self.properties_filter_description = ipw.HTML("<b>Filter by properties</b>:")

        self.job_state_dropdown = ipw.Dropdown(
            options={
                "Any": "",
                "Finished": "finished",
                "Running": "running",
                "Excepted": "except",
                "Killed": "killed",
            },
            value="",  # Default value corresponding to "Any"
            description="Job state:",
        )

        self.toggle_time_format = ipw.ToggleButtons(
            options=["Absolute", "Relative"],
            value="Absolute",  # Default to showing Absolute time
            description="Time format:",
        )
        self.toggle_time_format.observe(self.update_table_visibility, names="value")

        self.toggle_id_format = ipw.ToggleButtons(
            options=["pk", "uuid"],
            value="pk",
            description="ID format:",
        )
        self.toggle_id_format.observe(self.update_table_visibility, names="value")

        # Date pickers for range-based filtering
        self.time_start = ipw.DatePicker(description="Start time:")
        self.time_end = ipw.DatePicker(description="End time:")
        self.time_box = ipw.HBox([self.time_start, self.time_end])

        # Connect checkboxes and dropdowns to the filter logic
        for cb in property_checkboxes:
            cb.observe(self.apply_filters, names="value")
        self.time_start.observe(self.apply_filters, names="value")
        self.time_end.observe(self.apply_filters, names="value")
        self.job_state_dropdown.observe(self.apply_filters, names="value")

        display_options = ipw.VBox(
            children=[
                ipw.VBox(children=[self.toggle_time_format, self.toggle_id_format]),
            ],
            layout=ipw.Layout(
                border="1px solid lightgray",
                padding="0.5em",
            ),
        )

        filters = ipw.VBox(
            children=[
                ipw.VBox(
                    children=[
                        self.job_state_dropdown,
                        self.time_box,
                        ipw.VBox(
                            children=[
                                self.properties_filter_description,
                                self.properties_box,
                            ],
                            layout=ipw.Layout(
                                border="1px solid lightgray",
                                padding="0.5em",
                                margin="5px",
                            ),
                        ),
                    ]
                ),
            ],
            layout=ipw.Layout(
                border="1px solid lightgray",
                padding="0.5em",
            ),
        )

        self.main.children = [
            ipw.HTML("<h4>Display options:</h4>"),
            display_options,
            ipw.HTML("<h4>Filters:</h4>"),
            filters,
            self.table,
        ]

        self.update_table_value(self.data)

    def update_table_value(self, data_list):
        """Prepare the data to be shown (adding or hiding columns, etc.),
        and load it into `self.table`.
        """
        # Adjust which creation_time columns to hide
        COLUMNS["creation_time_relative"]["hide"] = (
            self.toggle_time_format.value != "Relative"
        )
        COLUMNS["creation_time_absolute"]["hide"] = (
            self.toggle_time_format.value != "Absolute"
        )

        # Build a new list that has an 'id' column, etc.
        display_data = []
        for row in data_list:
            row_copy = dict(row)
            # Switch the 'id' column depending on pk vs uuid toggle
            if self.toggle_id_format.value == "pk":
                row_copy["id"] = row_copy["pk_with_link"]
            else:
                row_copy["id"] = row_copy["uuid_with_link"]

            # Overwrite "state" with icon-based representation
            row_copy["state"] = determine_state_icon(row_copy)
            display_data.append(row_copy)

        # Figure out which columns to show the table
        columns = []
        for key in COLUMNS.keys():
            col_spec = dict(COLUMNS[key])
            col_spec["field"] = key
            columns.append(col_spec)

        self.table.from_data(display_data, columns=columns)

    def apply_filters(self, _):
        """Filter the raw data based on job state, selected properties,
        and date range. Then update the table display.
        """
        filtered = []

        selected_properties = [
            cb.description for cb in self.properties_box.children if cb.value
        ]

        # Convert DatePicker values (which are dates) into datetimes with UTC
        start_time = None
        end_time = None
        if self.time_start.value is not None:
            start_time = datetime.combine(self.time_start.value, datetime.min.time())
            start_time = start_time.replace(tzinfo=timezone.utc)

        if self.time_end.value is not None:
            end_time = datetime.combine(self.time_end.value, datetime.max.time())
            end_time = end_time.replace(tzinfo=timezone.utc)

        # State filter (empty string means "Any")
        desired_state_substring = self.job_state_dropdown.value

        for row in self.data:
            if desired_state_substring:
                if desired_state_substring not in row["state"].lower():
                    continue

            row_props = row["properties"]
            if selected_properties:
                # Must have all selected properties in row_props
                if not all(prop in row_props for prop in selected_properties):
                    continue

            ctime = row.get("creation_time", None)
            if ctime is not None:
                if start_time and ctime < start_time:
                    continue
                if end_time and ctime > end_time:
                    continue

            filtered.append(row)

        self.update_table_value(filtered)

    def update_table_visibility(self, _):
        """Called when toggles for time format or ID format change."""
        # simply re-apply filters (which triggers a re-draw).
        self.apply_filters(None)
