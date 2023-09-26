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
    # new workflow, this will reset the GUI
    app.work_chain_selector.value = None
    # check if the value are reload correctly
    assert app.configure_step.workchain_settings.relax_type.value == "positions_cell"
    assert app.configure_step.workchain_settings.spin_type.value == "none"
    assert app.configure_step.workchain_settings.properties["bands"].run.value is False
    assert app.configure_step.workchain_settings.properties["pdos"].run.value is False
