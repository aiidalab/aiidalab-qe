import pytest


@pytest.mark.usefixtures("sssp")
def test_settings(submit_app_generator):
    """Test the settings of the xas app."""
    app = submit_app_generator(properties=["xas"])
    configure_step = app.configure_step
    # test get_panel_value
    # select the first elmement
    configure_step.settings["xas"].element_and_ch_treatment.children[0].children[
        0
    ].value = True
    configure_step.settings["xas"].supercell_min_parameter.value = 4.0
    parameters = configure_step.settings["xas"].get_panel_value()
    assert parameters["core_hole_treatments"] == {"Si": "full"}
    assert parameters["pseudo_labels"] == {
        "Si": {
            "gipaw": "Si.pbe-van_gipaw.UPF",
            "core_hole": "Si.star1s-pbe-van_gipaw.UPF",
        }
    }
    assert parameters["core_wfc_data_labels"] == {"Si": "Si.pbe-van_gipaw.dat"}
    assert parameters["supercell_min_parameter"] == 4.0
    # test set_panel_value
    # update the parameters
    parameters["supercell_min_parameter"] = 5.0
    parameters["core_hole_treatments"] = {"Si": "xch_smear"}
    configure_step.settings["xas"].set_panel_value(parameters)
    assert configure_step.settings["xas"].supercell_min_parameter.value == 5.0
    assert (
        configure_step.settings["xas"]
        .element_and_ch_treatment.children[0]
        .children[1]
        .value
        == "xch_smear"
    )
