import ipywidgets as ipw
import pandas as pd
from IPython.display import display

from aiida.orm import QueryBuilder

state_icons = {
    "running": "⏳",
    "finished": "✅",
    "excepted": "⚠️",
    "killed": "❌",
}


class QueryInterface:
    def __init__(self):
        pass

    def setup_table(self):
        self.df = self.load_data()
        self.table = ipw.HTML()
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
            df["Creation time"] = df["ctime"].apply(
                lambda x: x.strftime("%Y-%m-%d %H:%M:%S")
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
                "Creation time",
                "Structure",
                "State",
                "Label",
                "Description",
                "Relax_type",
                "Delete",
                "Download",
                "Properties",
                "ctime",
            ]
        ]

    def setup_widgets(self):
        self.css_style = """
            <style>
                .df { border: none; }
                .df tbody tr:nth-child(odd) { background-color: #e5e7e9; }
                .df tbody tr:nth-child(odd):hover { background-color:   #f5b7b1; }
                .df tbody tr:nth-child(even):hover { background-color:  #f5b7b1; }
                .df tbody td { min-width: 80px; text-align: center; border: none }
                .df th { text-align: center; border: none;  border-bottom: 1px solid black;}
            </style>
            """

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
            children=property_checkboxes, description="Properties:"
        )
        self.properties_filter_description = ipw.HTML(
            "<p><b>Properties filter:</b> Select one or more properties to narrow the results. Only calculations that include all the selected properties will be displayed. Leave all checkboxes unselected to include calculations regardless of their properties.</p>"
        )
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
        self.label_search_field = ipw.Text(
            value="",
            placeholder="Enter a keyword",
            description="",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.label_search_description = ipw.HTML(
            "<p><b>Search label:</b> Enter a keyword to search in both the <i>Label</i> and <i>Description</i> fields. Matches will include any calculations where the keyword is found in either field.</p>"
        )
        self.toggle_description_checkbox = ipw.Checkbox(
            value=False,  # Show the Description column by default
            description="Show description",
            indent=False,
        )
        self.toggle_description_checkbox.observe(
            self.update_table_visibility, names="value"
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
        self.label_search_field.observe(self.apply_filters, names="value")

        self.filters_layout = ipw.VBox(
            [
                ipw.HTML("<h3>Search & filter calculations:</h3>"),
                ipw.VBox(
                    [
                        self.job_state_dropdown,
                        self.time_box,
                        ipw.HBox(
                            [
                                self.label_search_description,
                                self.label_search_field,
                            ]
                        ),
                        ipw.VBox(
                            [self.properties_filter_description, self.properties_box]
                        ),
                        #   self.apply_filters_btn,
                    ]
                ),
                ipw.HTML("<h3>Display options:</h3>"),
                ipw.VBox(
                    [
                        self.toggle_description_checkbox,
                        self.toggle_time_format,
                        self.toggle_id_format,
                    ]
                ),
            ]
        )
        self.get_table_value(self.df)

    def get_table_value(self, display_df):
        if display_df.empty:
            self.table.value = "<h2>No results found</h2>"
            return
        # Adjust the Creation time column based on the toggle state
        if self.toggle_time_format.value == "Relative":
            now = pd.Timestamp.now(tz="UTC")
            display_df["Creation time"] = display_df["ctime"].apply(
                lambda x: f"{(now - x).days}D ago" if pd.notnull(x) else "N/A"
            )
        else:
            display_df["Creation time"] = display_df["ctime"].apply(
                lambda x: x.strftime("%Y-%m-%d %H:%M:%S") if pd.notnull(x) else "N/A"
            )
        # Conditionally drop the Description column based on the checkbox state
        if not self.toggle_description_checkbox.value:
            display_df = display_df.drop(columns=["Description"])

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
            lambda x: f"{x.capitalize()}{state_icons.get(x.lower())}"
        )
        display_df.rename(
            columns={
                "State": "State 🟢",
                "Creation time": "Creation time ⏰",
                "ID": "ID 🔗",
                "Relax_type": "Relax type",
            },
            inplace=True,
        )
        display_df = display_df.drop(columns=["Properties", "ctime"])
        self.table.value = self.css_style + display_df.to_html(
            classes="df", escape=False, index=False
        )

    def apply_filters(self, _):
        selected_properties = [
            cb.description for cb in self.properties_box.children if cb.value
        ]
        filtered_df = self.df.copy()
        filtered_df = filtered_df[
            filtered_df["State"].str.contains(self.job_state_dropdown.value)
        ]
        if self.label_search_field.value:
            filtered_df = filtered_df[
                filtered_df["Label"].str.contains(
                    self.label_search_field.value, case=False, na=False
                )
                | filtered_df["Description"].str.contains(
                    self.label_search_field.value, case=False, na=False
                )
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
        self.get_table_value(filtered_df)

    def update_table_visibility(self, _):
        # Reapply filters to refresh the table visibility when the checkbox changes
        self.apply_filters(None)

    def display(self):
        display(self.filters_layout)
        display(self.table)
