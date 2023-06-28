def test_app_steps(app):
    """Test the link"""
    assert len(app.steps.steps) == 4
    #
    step1 = app.steps.steps[0][1]
    structure = step1.manager.children[0].children[3]
    structure.children[0].value = structure.children[0].options[1][1]
    step1.confirm()
    assert step1.state == step1.State.SUCCESS
    step2 = app.steps.steps[1][1]
    assert step2.previous_step_state == step1.State.SUCCESS


def test_app_submit(app):
    # Step 1: select structure from example
    step1 = app.steps.steps[0][1]
    structure = step1.manager.children[0].children[3]
    structure.children[0].value = structure.children[0].options[1][1]
    step1.confirm()
    # step 2
    step2 = app.steps.steps[1][1]
    step2.workchain_settings.relax_type.value = "none"
    step2.workchain_settings.properties["bands"].run.value = True
    step2.workchain_settings.properties["pdos"].run.value = True
    step2.basic_settings.workchain_protocol.value = "fast"
    step2.confirm()
    # step 3
    #
    step3 = app.steps.steps[2][1]
    assert step3.previous_step_state == step1.State.SUCCESS
    assert step3.state == step3.State.READY
    # max cpu is 2 in github runner
    step3.resources_config.num_cpus.value = 2
    assert step3.resources_config.num_cpus.value == 2
    builder, ui_parameters = step3._create_builder()
    # step3.submit()
