import pytest

from aiidalab_qe.app.main import App


@pytest.mark.usefixtures("aiida_profile_clean", "sssp")
def test_settings(submit_app_generator):
    """Test the settings of the xas app."""
    app: App = submit_app_generator(properties=["xas"])

    xas_model = app.configure_model.get_model("xas")
    xas_model.update()

    # Test getting the model state
    xas_model.elements["Si"] = True
    xas_model.supercell_min_parameter = 4.0
    parameters = xas_model.get_model_state()
    assert parameters["core_hole_treatments"] == {"Si": "full"}
    assert parameters["pseudo_labels"] == {
        "Si": {
            "gipaw": "Si.pbe-van_gipaw.UPF",
            "core_hole": "Si.star1s-pbe-van_gipaw.UPF",
        }
    }
    assert parameters["core_wfc_data_labels"] == {"Si": "Si.pbe-van_gipaw.dat"}
    assert parameters["supercell_min_parameter"] == 4.0

    # Test setting the model state
    parameters["supercell_min_parameter"] = 5.0
    parameters["core_hole_treatments"] = {"Si": "xch_smear"}
    xas_model.set_model_state(parameters)
    assert xas_model.supercell_min_parameter == 5.0
    assert xas_model.core_hole_treatments["Si"] == "xch_smear"
