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
checkpoint. The AiiDA daemon must be running for this action to succeed. When the
daemon is stopped, a warning appears at the top of the Results step and the Play
buttons are disabled until the daemon becomes active. When multiple paused
processes are listed, **Play all** resumes only the paused processes belonging to
the current workflow. The button is hidden when there are no paused processes.
Any error reported while requesting the resume is shown below the paused-process
table.

The reason text is supplied by AiiDA and can vary depending on the workflow and the
plugin that paused the process.

Terminating the workflow
=========================

If you do not want to resume a paused process, use the **Kill workflow** button at
the top of the Results step. This terminates the entire workflow rather than only
the selected paused process. If paused processes are present, the app first asks
AiiDA to resume them and waits for them to leave the paused state before sending
the kill request. The daemon must be running for this sequence to proceed.

Configuring failure behavior
=============================

The **On unhandled failure** setting is available under **Advanced settings**.
Choose one of the following behaviors:

* **Abort**: stop the workflow when an unhandled failure occurs.
* **Pause**: pause the workflow for inspection. This is the default.
* **Restart once**: retry once before aborting if the retry also fails.
* **Restart and pause**: retry once before pausing if the retry also fails.

The selected value is passed to the root app workflow and compatible first-level
plugin workflows. Plugin workflows are responsible for propagating the setting
to their own subworkflows.
