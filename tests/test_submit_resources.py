def test_resources_widget():
    """Test resources widget."""
    from aiidalab_qe.app.submit.resources import (
        ParallelizationSettings,
        ResourceSelectionWidget,
    )

    widget = ResourceSelectionWidget()
    widget.num_nodes.value = 2
    widget.num_cpus.value = 4
    assert widget.num_nodes.value == 2
    assert widget.num_cpus.value == 4
    widget.reset()
    assert widget.num_nodes.value == 1
    assert widget.num_cpus.value == 1
    #
    widget = ParallelizationSettings()
    widget.npools.value = 2
    assert widget.npools.value == 2
    widget.reset()
    assert widget.npools.value == 1


def test_resources():
    """Test get and set resources."""
    from aiidalab_qe.app.submit.submit import SubmitQeAppWorkChainStep

    step = SubmitQeAppWorkChainStep(qe_auto_setup=False)
    step.resources_config.num_nodes.value = 2
    step.resources_config.num_cpus.value = 4
    step.parallelization.npools.value = 2
    assert step.resources_config.num_nodes.value == 2
    assert step.resources_config.num_cpus.value == 4
    assert step.parallelization.npools.value == 2
    # test `get_resource` method
    resources = step.get_resource()
    assert resources["num_machines"] == 2
    assert resources["num_mpiprocs_per_machine"] == 4
    assert resources["npools"] == 2
