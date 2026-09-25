==================
Paused processes
==================

AiiDA processes can be paused when a workflow needs user intervention, for example
when a restart workflow reaches a configured iteration limit. A paused process keeps
its state and can be resumed from the app instead of restarting the workflow from
scratch.

Recognizing paused processes
============================

When one or more processes in the current workflow are paused, the Results step shows
a warning banner with the number of paused processes. The banner also points to the
Paused processes section in the Status view.

Reviewing a paused process
==========================

To inspect paused processes:

#. Open the **Status** view in the Results step.
#. Open **Paused Processes**.
#. Review the process primary key, label, and reason shown in the table.

The **Go to** action switches to the **Advanced View** and selects the corresponding
process in the workflow tree.

Resuming a process
==================

Use the **Play** action to ask AiiDA to resume the selected process from its saved
checkpoint. The AiiDA daemon must be running for this action to succeed. Any error
reported while requesting the resume is shown below the paused-process table.

The reason text is supplied by AiiDA and can vary depending on the workflow and the
plugin that paused the process.

Terminating the workflow
=========================

If you do not want to resume a paused process, use the **Kill workflow** button at
the top of the Results step. This terminates the entire workflow rather than only
the selected paused process.
