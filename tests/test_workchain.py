from aiida import engine, orm
from aiida.common import LinkType


def test_qeapp_workchain(generate_qeapp_workchain):
    wkchain = generate_qeapp_workchain()
    print("wkchain: ", wkchain)
    assert wkchain.setup() is None
    # generate structure for scf calculation
    wkchain.setup()
    #
    assert wkchain.should_run_relax() is True
    assert wkchain.should_run_bands() is True
    assert wkchain.should_run_pdos() is True
    # run relax and return the process
    relax_process = wkchain.run_relax()["workchain_relax"]
    # add test result
    result = orm.Dict({"energy": 0})
    result.store()
    result.base.links.add_incoming(
        relax_process, link_type=LinkType.RETURN, link_label="output_parameters"
    )
    # set state to finished.
    relax_process.set_process_state(engine.ProcessState.FINISHED)
    relax_process.set_exit_status(0)
    wkchain.ctx["workchain_relax"] = relax_process
