"""Regression tests for summaries of imported calculations."""

from uuid import uuid4

import pytest

from aiida import orm
from aiida.common.links import LinkType
from aiida.orm.utils.serialize import serialize
from aiidalab_qe.app.result.components.summary import WorkflowSummaryModel


@pytest.mark.parametrize(
    "pseudo_family",
    ["SSSP/1.3/PBE/efficiency", "PseudoDojo/0.4/PBE/SR/stringent/upf", None],
)
@pytest.mark.parametrize("missing", [False, True])
def test_summary_with_missing_pseudo(
    generate_structure_data, generate_upf_data, pseudo_family, missing
):
    """Missing UUIDs in saved UI settings must not prevent results rendering."""
    if missing:
        pseudo_uuid = str(uuid4())
    else:
        pseudo = generate_upf_data("Si").store()
        pseudo.base.extras.set_many({"functional": "PBE", "relativistic": "SR"})
        pseudo_uuid = str(pseudo.uuid)

    structure = generate_structure_data()
    process = orm.WorkChainNode()
    process.base.links.add_incoming(
        structure, link_type=LinkType.INPUT_WORK, link_label="structure"
    )
    process.store()
    parameters = {
        "workchain": {
            "relax_type": "none",
            "protocol": "moderate",
            "spin_type": "none",
            "electronic_type": "insulator",
        },
        "advanced": {
            "pseudo_family": pseudo_family,
            "kpoints_distance": 0.2,
            "pw": {
                "pseudos": {"Si": pseudo_uuid},
                "parameters": {"SYSTEM": {"ecutwfc": 30.0, "ecutrho": 240.0}},
            },
        },
    }
    saved_settings = serialize(parameters)
    process.base.extras.set("ui_parameters", saved_settings)
    model = WorkflowSummaryModel(process_uuid=process.uuid)

    report = model._generate_report_parameters()
    html = model.generate_report_html()

    assert report["advanced_settings"]["energy_cutoff_wfc"] == "30.0 Ry"
    if missing:
        assert "Unavailable" in report["advanced_settings"]["pseudos"][0]
        assert pseudo_uuid in html
        expected_functional = "PBE" if pseudo_family else "Unavailable"
    else:
        assert report["advanced_settings"]["pseudos"] == ["<b>Si:</b> Si.upf"]
        expected_functional = "PBE"
    assert report["advanced_settings"]["functional"]["value"] == expected_functional
    assert "30.0 Ry" in html
    assert process.base.extras.get("ui_parameters") == saved_settings
