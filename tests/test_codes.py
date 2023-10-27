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


def test_update_codes_display():
    """Test update_codes_display method.
    If the workchain property is not selected, the related code should be hidden.
    """
    from aiidalab_qe.app.submission import SubmitQeAppWorkChainStep

    submit = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    submit.update_codes_display()
    assert submit.codes["dos"].layout.display == "none"
    submit.input_parameters = {"workchain": {"properties": ["pdos"]}}
    submit.update_codes_display()
    assert submit.codes["dos"].layout.display == "block"


def test_identify_submission_blockers(app):
    """Test identify_submission_blockers method."""
    submit = app.submit_step
    blockers = list(submit._identify_submission_blockers())
    # there is one blocker: ['The SSSP library is not installed.']
    assert len(blockers) == 1
    submit.input_parameters = {"workchain": {"properties": ["pdos"]}}
    blockers = list(submit._identify_submission_blockers())
    assert len(blockers) == 1
    # set dos code to None, will introduce another blocker
    dos_value = submit.codes["dos"].value
    submit.codes["dos"].value = None
    blockers = list(submit._identify_submission_blockers())
    assert len(blockers) == 2
    # set dos code back will remove the blocker
    submit.codes["dos"].value = dos_value
    blockers = list(submit._identify_submission_blockers())
    assert len(blockers) == 1
