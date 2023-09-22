def test_reload(submit_app_generator, generate_qeapp_workchain):
    """Test if the GUI paramters can be reload properly"""
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
