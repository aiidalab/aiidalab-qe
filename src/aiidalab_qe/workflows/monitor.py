import tempfile

from aiida.orm import CalcJobNode
from aiida.transports import Transport


def relax_pw(node: CalcJobNode, transport: Transport) -> str:
    """Retrieve and inspect files in working directory of job to determine whether the job should be killed.

    :param node: The node representing the calculation job.
    :param transport: The transport that can be used to retrieve files from remote working directory.
    :returns: A string if the job should be killed, `None` otherwise.
    """
    if "relax" not in node.base.links.get_outgoing().all_link_labels():
        return
    node = node.base.links.get_outgoing().get_node_by_label("relax")
    with tempfile.NamedTemporaryFile("w+") as handle:
        transport.getfile(node.options.output_filename, handle.name)
        handle.seek(0)
        output = handle.read()

    energies = []
    for line in output:
        if "total energy              =" in line:
            energies.append(float(line.split("=")[1]))
    node.base.extras.set("monitor.energies", energies)
