def test_code_not_selected(submit_app_generator):
    """Test if there is an error when the code is not selected."""
    app = submit_app_generator()
    app.submit_step.codes["dos"].value = None
    app.submit_step._create_builder()


def test_set_selected_codes(submit_app_generator):
    """Test set_selected_codes method."""
    from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep

    app = submit_app_generator()
    submit_step = app.submit_step

    submit_step._create_builder()

    new_submit_step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    new_submit_step.set_selected_codes(submit_step.ui_parameters["codes"])

    assert new_submit_step.get_selected_codes() == submit_step.get_selected_codes()


def test_set_code_status():
    """Test set_codes_status method.
    If the workchain property is not selected, the related code should be disabled.
    """
    from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep

    submit = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    submit.set_codes_status()
    assert submit.codes["dos"].code_select_dropdown.disabled is True
    submit.input_parameters = {"workchain": {"properties": ["pdos"]}}
    submit.set_codes_status()
    assert submit.codes["dos"].code_select_dropdown.disabled is False
