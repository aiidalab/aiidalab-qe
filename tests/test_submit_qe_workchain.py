import pytest


@pytest.mark.usefixtures("aiida_profile_clean")
def test_get_input_parameters(data_regression):
    from aiidalab_qe.pseudos import PseudoFamilySelector
    from aiidalab_qe.steps import (
        KpointSettings,
        SmearingSettings,
        SubmitQeAppWorkChainStep,
        WorkChainSettings,
    )

    submit_step = SubmitQeAppWorkChainStep(qe_auto_setup=False)

    submit_step.workchain_settings = WorkChainSettings()
    submit_step.pseudo_family_selector = PseudoFamilySelector()
    submit_step.kpoints_settings = KpointSettings()
    submit_step.smearing_settings = SmearingSettings()

    parameters = submit_step.get_input_parameters()

    # Check and validate the parameters
    data_regression.check(parameters)
