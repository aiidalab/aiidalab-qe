def test_xps(app):
    step1 = app.structure_step
    structure = step1.manager.children[0].children[3]
    # select water molecule
    structure.children[0].value = structure.children[0].options[-1][1]
    step1.confirm()
    # Step 2: configure calculation
    step2 = app.configure_step
    step2.workchain_settings.properties["xps"].run.value = True
    # test the peak list are udpate
    assert len(step2.settings["xps"].peak_list.children) == 2
    # selet the peak
    step2.settings["xps"].peak_list.children[0].value = True
    step2.confirm()
