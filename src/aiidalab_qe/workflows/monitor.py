import tempfile

from aiida.orm import CalcJobNode
from aiida.transports import Transport


def relax_pw(node: CalcJobNode, transport: Transport) -> str:
    """Retrieve and save the total energy from aiida.out file and save it into the node extras.

    :param node: The node representing the calculation job.
    :param transport: The transport that can be used to retrieve files from remote working directory.
    :returns: A string if the job should be killed, `None` otherwise.
    """
    try:
        with tempfile.NamedTemporaryFile("w+") as handle:
            transport.getfile(node.get_option("output_filename"), handle.name)
            handle.seek(0)
            output = handle.read()

        energies = []
        lines = output.splitlines()
        for line in lines:
            if "total energy              =" in line:
                energies.append(float(line.split("=")[1].split()[0]))
        node.base.extras.set("monitor_energies", energies)
    except FileExistsError:
        return
