def test_reload_and_reset(submit_app_generator, generate_qeapp_workchain):
    """Test if the GUI paramters can be reload and reset properly"""
    wkchain = generate_qeapp_workchain(
        relax_type="positions", run_bands=True, run_pdos=False, spin_type="collinear"
    )
    app = submit_app_generator()
    # select the pk
    app.work_chain_selector.value = wkchain.node.pk
    # check if the value are reload correctly
    assert app.configure_step.workchain_settings.relax_type.value == "positions"
    assert app.configure_step.workchain_settings.spin_type.value == "collinear"
    assert app.configure_step.workchain_settings.properties["bands"].run.value is True
    assert app.configure_step.workchain_settings.properties["pdos"].run.value is False
    assert (
        len(
            app.configure_step.advanced_settings.pseudo_setter.pseudo_setting_widgets.children
        )
        > 0
    )
    # new workflow, this will reset the GUI
    app.work_chain_selector.value = None
    # check if the value are reload correctly
    assert app.structure_step.manager.structure is None
    assert app.configure_step.workchain_settings.relax_type.value == "positions_cell"
    assert app.configure_step.workchain_settings.spin_type.value == "none"
    assert app.configure_step.workchain_settings.properties["bands"].run.value is False
    assert app.configure_step.workchain_settings.properties["pdos"].run.value is False
    assert app.configure_step.advanced_settings.pseudo_setter.ecutwfc_setter.value == 0
    assert (
        len(
            app.configure_step.advanced_settings.pseudo_setter.pseudo_setting_widgets.children
        )
        == 0
    )
    assert app.submit_step.resources_config.num_cpus.value == 1


def test_select_new_structure(app_to_submit, generate_structure_data):
    """Test if the new structure is selected, the confirmed structure is reset"""
    app = app_to_submit
    assert app.structure_step.confirmed_structure is not None
    # select a new structure will reset the confirmed structure
    app.structure_step.structure = generate_structure_data()
    assert app.structure_step.confirmed_structure is None


def test_unsaved_changes(app_to_submit):
    """Test if the unsaved changes are handled correctly"""
    from aiidalab_widgets_base import WizardAppWidgetStep

    app = app_to_submit
    # go to the configue step, and make some changes
    app._wizard_app_widget.selected_index = 1
    app.configure_step.workchain_settings.relax_type.value = "positions"
    # go to the submit step
    app._wizard_app_widget.selected_index = 2
    # the state of the configue step should be updated.
    assert app.configure_step.state == WizardAppWidgetStep.State.CONFIGURED
    # check if a new blocker is added
    assert len(app.submit_step.external_submission_blockers) == 1
    # confirm the changes
    app._wizard_app_widget.selected_index = 1
    app.configure_step.confirm()
    app._wizard_app_widget.selected_index = 2
    # the blocker should be removed
    assert len(app.submit_step.external_submission_blockers) == 0
