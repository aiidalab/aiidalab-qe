import ipywidgets as ipw
import pandas as pd
from table_widget import TableWidget

from aiida.orm import QueryBuilder, load_node
from aiidalab_qe.common.widgets import LoadingWidget

STATE_ICONS = {
    "running": "⏳",
    "finished": "✅",
    "excepted": "⚠️",
    "killed": "❌",
}

COLUMNS = {
    "ID": {"headerName": "ID 🔗", "dataType": "link", "editable": False},
    "Creation time absolute": {
        "headerName": "Creation Time ⏰ (absolute)",
        "type": "date",
        "width": 150,
        "editable": False,
    },
    "Creation time relative": {
        "headerName": "Creation Time ⏰ (relative)",
        "width": 150,
        "editable": False,
    },
    "Structure": {"headerName": "Structure", "editable": False},
    "State": {"headerName": "State 🟢", "editable": False},
    "Label": {"headerName": "Label", "width": 300, "editable": True},
    "Description": {
        "headerName": "Description",
        "width": 300,
        "editable": True,
        "hide": True,
    },
    "Relax_type": {"headerName": "Relax_type", "editable": False, "hide": True},
    "Delete": {"headerName": "Delete", "dataType": "link", "editable": False},
    "Download": {"headerName": "Download", "dataType": "link", "editable": False},
    "UUID": {"headerName": "UUID", "editable": False, "hide": True},
}


class CalculationHistory:
    def __init__(self):
        self.main = ipw.VBox(children=[LoadingWidget("Loading the table...")])

        self.table = TableWidget(style={"margin-top": "20px"})

        def on_row_update(change):
            node = load_node(change["new"]["UUID"])
            node.label = change["new"]["Label"]
            node.description = change["new"]["Description"]

        self.table.observe(on_row_update, "updatedRow")

    def load_table(self):
        self.df = self.load_data()
        self.setup_widgets()

    def load_data(self):
        from aiidalab_qe.workflows import QeAppWorkChain

        projections = [
            "id",
            "uuid",
            "extras.structure",
            "ctime",
            "attributes.process_state",
            "label",
            "description",
            "extras.workchain.relax_type",
            "extras.workchain.properties",
        ]
        headers = [
            "PK",
            "UUID",
            "Structure",
            "ctime",
            "State",
            "Label",
            "Description",
            "Relax_type",
            "Properties",
        ]

        qb = QueryBuilder()
        qb.append(QeAppWorkChain, project=projections, tag="process")
        qb.order_by({"process": {"ctime": "desc"}})
        results = qb.all()

        df = pd.DataFrame(results, columns=headers)
        # Check if DataFrame is not empty
        if not df.empty:
            df["Creation time absolute"] = df["ctime"].apply(
                lambda x: x.strftime("%Y-%m-%d %H:%M:%S")
            )
            now = pd.Timestamp.now(tz="UTC")
            df["Creation time relative"] = df["ctime"].apply(
                lambda x: f"{(now - x).days}D ago" if pd.notnull(x) else "N/A"
            )
            df["Delete"] = df["PK"].apply(
                lambda pk: f'<a href="./delete.ipynb?pk={pk}" target="_blank">Delete</a>'
            )
            df["Download"] = df["PK"].apply(
                lambda pk: f'<a href="./download.ipynb?pk={pk}" target="_blank">Download</a>'
            )
            # add a link to the pk so that the user can inspect the calculation
            df["PK_with_link"] = df["PK"].apply(
                lambda pk: f'<a href="./qe.ipynb?pk={pk}" target="_blank">{pk}</a>'
            )
            # Store initial part of the UUID
            df["UUID_with_link"] = df.apply(
                lambda row: f'<a href="./qe.ipynb?pk={row["PK"]}" target="_blank">{row["UUID"][:8]}</a>',
                axis=1,
            )
            # replace all "waiting" states with "running"
            df["State"] = df["State"].apply(
                lambda x: "running" if x == "waiting" else x
            )
        else:
            # Initialize empty columns for an empty DataFrame
            df["Creation time"] = pd.Series(dtype="str")
            df["Delete"] = pd.Series(dtype="str")
        return df[
            [
                "PK_with_link",
                "UUID_with_link",
                "Creation time absolute",
                "Creation time relative",
                "Structure",
                "State",
                "Label",
                "Description",
                "Relax_type",
                "UUID",
                "Delete",
                "Download",
                "Properties",
                "ctime",
            ]
        ]

    def setup_widgets(self):
        unique_properties = set(self.df["Properties"].explode().dropna())
        unique_properties.discard(None)
        property_checkboxes = [
            ipw.Checkbox(
                value=False,
                description=prop,
                Layout=ipw.Layout(description_width="initial"),
                indent=False,
            )
            for prop in unique_properties
        ]
        self.properties_box = ipw.HBox(
            children=property_checkboxes,
            indent=True,
        )
        self.properties_filter_description = ipw.HTML("<b>Filter by properties</b>:")
        # Replace 'None' in 'Properties' with an empty list
        self.df["Properties"] = self.df["Properties"].apply(
            lambda x: [] if x is None else x
        )
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
            value="Absolute",  # Default to Absolute time
            description="Time format:",
        )
        self.toggle_time_format.observe(self.update_table_visibility, names="value")
        self.toggle_id_format = ipw.ToggleButtons(
            options=["PK", "UUID"],
            value="PK",  # Default to PK
            description="ID format:",
        )
        self.toggle_id_format.observe(self.update_table_visibility, names="value")
        # self.toggle_multi_selection = ipw.ToggleButtons(
        #     options=["On", "Off"],
        #     value="Off",
        #     description="Checkbox selection",
        #     tooltips=["Enable multiple selection.", "Disable multiple selection"],
        # )
        # self.toggle_multi_selection.observe(
        #     self.update_table_configuration, names="value"
        # )

        self.time_start = ipw.DatePicker(description="Start time:")
        self.time_end = ipw.DatePicker(description="End time:")
        self.time_box = ipw.HBox([self.time_start, self.time_end])
        # self.apply_filters_btn = ipw.Button(description='Apply Filters')
        # self.apply_filters_btn.on_click(self.apply_filters)
        for cb in property_checkboxes:
            cb.observe(self.apply_filters, names="value")
        self.time_start.observe(self.apply_filters, names="value")
        self.time_end.observe(self.apply_filters, names="value")
        self.job_state_dropdown.observe(self.apply_filters, names="value")
        display_options = ipw.VBox(
            children=[
                ipw.HTML("<h4>Display options:</h4>"),
                ipw.VBox(
                    [
                        self.toggle_time_format,
                        self.toggle_id_format,
                        # self.toggle_multi_selection,
                    ]
                ),
            ],
            layout=ipw.Layout(
                border="1px solid #ddd",
                padding="0px",
                margin="0px, 0px, 100px, 0px",
                background_color="#f9f9f9",
                border_radius="5px",
                box_shadow="0 2px 5px rgba(0, 0, 0, 0.1)",
            ),
        )
        filters = ipw.VBox(
            children=[
                ipw.HTML("<h4>Filters:</h4>"),
                ipw.VBox(
                    [
                        self.job_state_dropdown,
                        self.time_box,
                        ipw.VBox(
                            [self.properties_filter_description, self.properties_box],
                            layout=ipw.Layout(
                                border="1px solid #ddd", padding="5px", margin="5px"
                            ),
                        ),
                    ]
                ),
            ],
            layout=ipw.Layout(
                border="1px solid #ddd",
                padding="0px",
                margin="0px, 0px, 100px, 0px",
                background_color="#f9f9f9",
                border_radius="5px",
                box_shadow="0 2px 5px rgba(0, 0, 0, 0.1)",
            ),
        )

        self.main.children = [
            display_options,
            filters,
            self.table,
        ]
        self.update_table_value(self.df)

    def update_table_value(self, display_df):
        # Adjust the Creation time column based on the toggle state
        COLUMNS["Creation time relative"]["hide"] = (
            self.toggle_time_format.value != "Relative"
        )
        COLUMNS["Creation time absolute"]["hide"] = (
            self.toggle_time_format.value != "Absolute"
        )
        # Adjust the ID column based on the toggle state
        if self.toggle_id_format.value == "PK":
            display_df = display_df.rename(columns={"PK_with_link": "ID"}).drop(
                columns=["UUID_with_link"]
            )
        else:
            display_df = display_df.rename(columns={"UUID_with_link": "ID"}).drop(
                columns=["PK_with_link"]
            )

        #
        display_df["State"] = display_df["State"].apply(
            lambda x: f"{x.capitalize()}{STATE_ICONS.get(x.lower())}"
        )
        display_df = display_df.drop(columns=["Properties", "ctime"])
        columns = []
        for col in display_df.columns:
            column = COLUMNS.get(col)
            column["field"] = col
            columns.append(COLUMNS[col])
        self.table.from_data(display_df, columns=columns)

    def apply_filters(self, _):
        selected_properties = [
            cb.description for cb in self.properties_box.children if cb.value
        ]
        filtered_df = self.df.copy()
        filtered_df = filtered_df[
            filtered_df["State"].str.contains(self.job_state_dropdown.value)
        ]
        if selected_properties:
            filtered_df = filtered_df[
                filtered_df["Properties"].apply(
                    lambda x: all(item in x for item in selected_properties)
                )
            ]
        if self.time_start.value and self.time_end.value:
            start_time = pd.to_datetime(self.time_start.value).normalize()
            end_time = pd.to_datetime(self.time_end.value).normalize() + pd.Timedelta(
                days=1, milliseconds=-1
            )
            start_time = start_time.tz_localize("UTC")
            end_time = end_time.tz_localize("UTC")
            filtered_df = filtered_df[
                (filtered_df["ctime"] >= start_time)
                & (filtered_df["ctime"] <= end_time)
            ]
        self.update_table_value(filtered_df)

    def update_table_visibility(self, _):
        # Reapply filters to refresh the table visibility when the checkbox changes
        self.apply_filters(None)

    def update_table_configuration(self, _):
        enable = True if self.toggle_multi_selection.value == "On" else False
        config = {"pagination": True, "checkboxSelection": enable}
        self.table.config = config
