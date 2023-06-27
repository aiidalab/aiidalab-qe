def test_app(app):
    assert len(app.steps.steps) == 4


def test_app_submit(app):
    from aiidalab_widgets_base import WizardAppWidgetStep

    was = WizardAppWidgetStep()
    # Step 1: select structure from example
    s1 = app.steps.steps[0][1]
    structure = s1.manager.children[0].children[3]
    structure.children[0].value = structure.children[0].options[1][1]
    s1.confirm()
    # step 2
    s2 = app.steps.steps[1][1]
    s2.workchain_settings.relax_type.value = "none"
    s2.workchain_settings.properties["bands"].run.value = True
    s2.workchain_settings.properties["pdos"].run.value = True
    s2.basic_settings.workchain_protocol.value = "fast"
    parameters = s2.get_input_parameters()
    s2.confirm()
    print("parameters: ", parameters)
    # step 3
    #
    s3 = app.steps.steps[2][1]
    assert s3.previous_step_state == was.State.SUCCESS
    assert s3.state == s3.State.CONFIGURED
    s3.resources_config.num_cpus.value = 4
    assert s3.resources_config.num_cpus.value == 4
    builder = s3.get_builder()
    print(builder)
    s3.submit()
