import pytest
from bs4 import BeautifulSoup

from aiidalab_qe.app.result import ViewQeAppWorkChainStatusAndResultsStep
from aiidalab_qe.app.result.components.summary import WorkChainSummaryModel
from aiidalab_qe.app.result.components.viewer import (
    WorkChainResultsViewer,
    WorkChainResultsViewerModel,
)
from aiidalab_qe.app.result.components.viewer.structure import (
    StructureResultsModel,
    StructureResultsPanel,
)
from aiidalab_qe.app.wizard_app import WizardApp


def test_result_step(app_to_submit, generate_qeapp_workchain):
    """Test the result step is properly updated when the process
    is running."""
    app: WizardApp = app_to_submit
    step: ViewQeAppWorkChainStatusAndResultsStep = app.results_step
    model = app.results_model
    model.process_uuid = generate_qeapp_workchain().node.uuid
    assert step.state == step.State.ACTIVE
    step.render()
    assert step.toggle_controls.value == "Status"
    results_model: WorkChainResultsViewerModel = model.get_model("results")  # type: ignore
    # All jobs are completed, so there should be no process status notifications
    for _, model in results_model.get_models():
        assert model.process_status_notification == ""


def test_kill_and_clean_buttons(app_to_submit, generate_qeapp_workchain):
    """Test the kill and clean_scratch button are properly displayed when the process
    is in different states."""
    step = app_to_submit.results_step
    step.render()
    model = app_to_submit.results_model
    model.process_uuid = generate_qeapp_workchain().node.uuid
    assert step.kill_button.layout.display == "block"
    assert step.clean_scratch_button.layout.display == "none"


def test_workchainview(generate_qeapp_workchain):
    """Test the result tabs are properly updated"""
    workchain = generate_qeapp_workchain()
    workchain.node.seal()
    model = WorkChainResultsViewerModel()
    viewer = WorkChainResultsViewer(model=model)
    model.process_uuid = workchain.node.uuid
    viewer.render()
    assert len(viewer.tabs.children) == 2
    assert viewer.tabs._titles["0"] == "Structure"  # type: ignore


def test_summary_report(data_regression, generate_qeapp_workchain):
    """Test the summary report can be properly generated."""
    workchain = generate_qeapp_workchain()
    model = WorkChainSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_parameters = model._generate_report_parameters()
    # Discard variable parameters
    for key in (
        "pk",
        "uuid",
        "creation_time",
        "modification_time",
    ):
        report_parameters["workflow_properties"].pop(key)
    for key in (
        "structure_pk",
        "structure_uuid",
    ):
        report_parameters["initial_structure_properties"].pop(key)
    data_regression.check(report_parameters)


def test_summary_report_advanced_settings(data_regression, generate_qeapp_workchain):
    """Test advanced settings are properly reported"""
    workchain = generate_qeapp_workchain(
        spin_type="collinear", electronic_type="metal", initial_magnetic_moments=0.1
    )
    model = WorkChainSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_parameters = model._generate_report_parameters()
    moments = report_parameters["advanced_settings"]["initial_magnetic_moments"]
    assert moments["Si"] == 0.1


@pytest.mark.parametrize(
    ("pbc", "symmetry_key"),
    [
        [(False, False, False), "point_group"],  # 0D
        [(True, False, False), "space_group"],  # 1D
        [(True, True, False), "space_group"],  # 2D
        [(True, True, True), "space_group"],  # 3D
    ],
)
def test_summary_report_symmetry_group(
    generate_qeapp_workchain,
    generate_structure_data,
    pbc,
    symmetry_key,
):
    """Test summary report includes correct symmetry group for all system dimension."""

    system = generate_structure_data("silicon", pbc=pbc)
    workchain = generate_qeapp_workchain(
        structure=system,
        run_bands=False,
        relax_type="none",
    )
    model = WorkChainSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_parameters = model._generate_report_parameters()
    assert symmetry_key in report_parameters["initial_structure_properties"]


def test_summary_view(generate_qeapp_workchain):
    """Test the report html can be properly generated."""
    workchain = generate_qeapp_workchain()
    model = WorkChainSummaryModel()
    model.process_uuid = workchain.node.uuid
    report_html = model.generate_report_html()
    parsed = BeautifulSoup(report_html, "html.parser")
    parameters = {
        "Energy cutoff (wave functions)": "30.0 Ry",
        "Total charge": "0.0",
    }
    for key, value in parameters.items():
        key_td = parsed.find("td", string=lambda tag, key=key: tag and key in tag.text)
        value_td = key_td.find_next_sibling("td")
        assert value in value_td.text


def test_structure_results_panel(generate_qeapp_workchain):
    """Test the structure results panel can be properly generated."""

    model = StructureResultsModel()
    panel = StructureResultsPanel(model=model)

    def test_table_data(model):
        rows = model.table_data[1:]  # skip table header
        for i, row in enumerate(rows):
            position = model.structure.sites[i].position
            x, y, z = (f"{coordinate:.2f}" for coordinate in position)
            assert row == [i + 1, "Si", 0, x, y, z]  # type: ignore

    assert model.title == "Structure"

    wc = generate_qeapp_workchain(relax_type="none")
    model.process_uuid = wc.node.uuid
    node = model.fetch_process_node()
    assert "Si<sub>2</sub>" in model.header
    assert "Initial" in model.sub_header
    assert "properties" in model.source  # inputs
    assert model.structure.pk == node.inputs.structure.pk
    assert str(node.inputs.structure.pk) in model.info
    test_table_data(model)

    panel.render()
    assert panel.view_toggle_button.layout.display == "none"

    wc = generate_qeapp_workchain(relax_type="positions_cell")
    model.process_uuid = wc.node.uuid
    node = model.fetch_process_node()
    assert "Initial" in model.sub_header
    assert panel.view_toggle_button.layout.display == "block"
    assert panel.view_toggle_button.description == "View relaxed"
    panel.view_toggle_button.click()
    assert panel.view_toggle_button.description == "View initial"
    assert "Relaxed" in model.sub_header
    assert "properties" not in model.source  # outputs
    assert model.structure.pk == node.outputs.structure.pk
    assert str(node.outputs.structure.pk) in model.info
    test_table_data(model)
