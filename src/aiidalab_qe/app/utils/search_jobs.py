import ipywidgets as ipw
import pandas as pd
from aiida.orm import QueryBuilder
from aiidalab_qe.workflows import QeAppWorkChain
from IPython.display import display


class QueryInterface:
    def __init__(self):
        self.df = self.load_data()
        self.table = ipw.HTML()
        self.setup_widgets()

    def load_data(self):
        projections = [
            "id",
            "extras.structure",
            "ctime",
            "attributes.process_state",
            "label",
            "extras.workchain.relax_type",
            "extras.workchain.properties",
        ]
        headers = [
            "PK",
            "Structure",
            "ctime",
            "State",
            "Label",
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
            df["Inspect"] = df["PK"].apply(
                lambda pk: f'<a href="./qe.ipynb?pk={pk}" target="_blank">Inspect</a>'
            )
        else:
            # Initialize empty columns for an empty DataFrame
            df["Creation time"] = pd.Series(dtype="str")
            df["Delete"] = pd.Series(dtype="str")
            df["Inspect"] = pd.Series(dtype="str")
        return df[
            [
                "PK",
                "Creation time",
                "Structure",
                "State",
                "Label",
                "Relax_type",
                "Delete",
                "Inspect",
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

        unique_properties = set(self.df["Properties"].explode())
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
        # Replace 'None' in 'Properties' with an empty list
        self.df["Properties"] = self.df["Properties"].apply(
            lambda x: [] if x is None else x
        )
        self.job_state_dropdown = ipw.Dropdown(
            options=["", "finished", "waiting", "except", "killed"],
            value="",
            description="Job State:",
        )
        self.label_search_field = ipw.Text(
            value="",
            placeholder="Enter label to search",
            description="Search Label:",
            disabled=False,
            style={"description_width": "initial"},
        )
        self.time_start = ipw.DatePicker(description="Start Time:")
        self.time_end = ipw.DatePicker(description="End Time:")
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
                ipw.HTML("<h2>Search results:</h2>"),
                ipw.VBox(
                    [
                        ipw.HBox(
                            [ipw.HTML("<h>Properties: </h>"), self.properties_box]
                        ),
                        self.label_search_field,
                        self.job_state_dropdown,
                        self.time_box,
                        #   self.apply_filters_btn,
                    ]
                ),
            ]
        )
        self.get_table_value(self.df)

    def get_table_value(self, display_df):
        if display_df.empty:
            self.table.value = "<h2>No results found</h2>"
            return
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

    def display(self):
        display(self.filters_layout)
        display(self.table)
