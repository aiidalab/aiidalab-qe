from aiidalab_qe.common.widgets import TableWidget


def test_table_widget_loads_assets_and_syncs_traits():
    table = TableWidget(data=[["Element"], ["Si"]], selected_rows=[0])

    assert "export const render" in table._esm
    assert ".custom-table table" in table._css
    assert table.data == [["Element"], ["Si"]]
    assert table.selected_rows == [0]

    table.selected_rows = []
    assert table.selected_rows == []
