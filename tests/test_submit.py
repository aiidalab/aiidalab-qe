def test_codes(pw_code, dos_code, projwfc_code, pw_code_70):
    """Test get and set codes."""
    from aiida.orm import load_code

    from aiidalab_qe.app.parameters import DEFAULT_PARAMETERS
    from aiidalab_qe.app.submit import SubmitQeAppWorkChainStep

    step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    # get the default codes
    assert (
        load_code(step.pw_code.value).full_label
        == DEFAULT_PARAMETERS["codes"]["pw_code"]
    )
    assert (
        load_code(step.dos_code.value).full_label
        == DEFAULT_PARAMETERS["codes"]["dos_code"]
    )
    assert (
        load_code(step.projwfc_code.value).full_label
        == DEFAULT_PARAMETERS["codes"]["projwfc_code"]
    )
    # test `get_selected_codes` method
    codes = step.get_selected_codes()
    codes["pw_code"] = pw_code.uuid
    # set the codes manually
    step.pw_code.value = pw_code_70.uuid
    assert load_code(step.pw_code.value).full_label == pw_code_70.full_label
    # set `set_selected_codes` method
    step.set_selected_codes(parameters=DEFAULT_PARAMETERS)
    assert step.pw_code.value == pw_code.uuid


def test_create_builder_default(
    app_to_submit,
    data_regression,
):
    """ "Test the creation of the workchain builder.

    metal, non-magnetic
    """

    submit_step = app_to_submit.steps.steps[2][1]

    builder, ui_parameters = submit_step._create_builder()

    # check and validate the builder
    got = builder_to_readable_dict(builder._inputs())

    # regression test
    data_regression.check(got)


def builder_to_readable_dict(builder):
    """transverse the builder and return a dictionary with readable values."""
    from aiida import orm
    from aiida.engine import ProcessBuilderNamespace

    ignore_keys = ["metadata", "monitors", "pseudos", "code", "structure"]

    readable_dict = {}
    for k, v in builder.items():
        if k in ignore_keys:
            continue
        elif isinstance(v, (dict, ProcessBuilderNamespace)):
            readable_dict[k] = builder_to_readable_dict(v)
        elif isinstance(v, orm.Dict):
            readable_dict[k] = v.get_dict()
        elif isinstance(v, (orm.Int, orm.Float, orm.Str, orm.Bool)):
            readable_dict[k] = v.value
        elif isinstance(v, orm.List):
            readable_dict[k] = v.get_list()
        else:
            readable_dict[k] = v

    return readable_dict
