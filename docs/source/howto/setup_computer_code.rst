===================================
How to setup a remote computer/code
===================================

The computational resources required for a calculation will often surpass those available on a local machine. In such cases, a user may want to leverage the computing power of a supercomputing cluster. The following guide will show you how to set up the AiiDAlab environment to run calculations on such remote machines.

.. note::

   You may notice on this page dropdown sections pertaining to "**CSCS**"

   .. dropdown:: ... CSCS

      A note regarding CSCS...

   These collapsible sections provide instructions specific to machines housed at the Swiss Supercomputer Centre (`CSCS`_), such as `Daint` and `Eiger`. However, if you do not have an account with CSCS, you may still find these instructions useful as an example for your machine.

----

Please expand **Step 3: Choose computational resources**.

.. figure:: /_static/images/step3.png
   :width: 12cm

Under the `Codes` section, you can set up executables for the three relevant Quantum ESPRESSO codes used in the app:

*  pw.x -  used in planewave (PW) self-consistent calculations
*  dos.x - used to compute the density-of-states (DOS)
*  projwfc.x - used to compute projected density-of-states (PDOS)

.. note::

   Code selection for density-of-states (DOS/PDOS) will be disabled unless **Projected density of states** is selected as a property of instance in step 2.

Codes are represented as ``code-label@computer-label``, where a `computer` refers to a collection of configurations relating to the remote machine. So we first need to define a computer (on the remote machine), then a code (an executable) on the remote machine attached to that computer.

To begin, click on `Setup new code` for pw.x. This will open the setup widget explained in the following section.

Setup Widget
------------

The Quantum ESPRESSO app includes a built-in widget meant to assist with setting up new codes/computers. We will show you how to use it here to set up the AiiDAlab environment with executable codes located on a remote machine.

.. figure:: /_static/images/setup_widget.png
   :width: 12cm

.. important::

   If you're using a CSCS machine, please check out the dropdown below. Otherwise, click on the `Detailed Setup` tab and proceed to the next section.

.. dropdown:: Predefined CSCS Machines

   The Quantum ESPRESSO app comes pre-configured with default codes located on various CSCS machines. We provide here a less-detailed, bulletpoint-styled guide to set up the environment for submitting calculations to a CSCS remote machine. For more details on each step, see the sections below.

   .. warning::

      If the dropdown list in step 2 below does not yet include a computer associated with your account/project (e.g. multicore-em01 for project em01), please click on the `Detailed Setup` tab and proceed with the sections below to establish a connection to the remote machine, and set up and define codes on a computer associated with it.

   **Steps:**

   #. Select a domain (e.g. daint.cscs.ch)
   #. Select a computer (e.g. multicore-em01)
   #. Select a code (e.g. pw-7.0-multicore)
   #. In `Quick Setup`, enter your username and click `Quick Setup`
   #. Enter your password in the **passwordless enabling log** box and press `Continue`
   #. Repeat steps 1-3 for each of the codes required for your calculation

   Once a code is created, you may select it from its respective dropdown menu.

Connection
^^^^^^^^^^

We need to set up a connection to the remote host. In **Set up password-less SSH connection**, you can set up an SSH connection including a proxy jump and authentication options specific to your supercomputing center of choice.

.. figure:: /_static/images/passwordless.png
   :width: 12cm

.. important::

   You will only need to do this once per host, i.e. computers set up on this host in the following section will all have access to the connection defined here.

.. note::

   If you are unsure of the details of the connection, you may consider reaching out to your institution's IT department or supercomputing IT service desk.

Once you've entered all the necessary details, click `Setup ssh`.

.. tip::

   You can test your connection by entering your password in the **passwordless enabling log** box above and clicking `Continue`. After a few second, you'll be asked to reenter your password. Please do so and press `Continue` again. If you see the message ``30s timeout``, please click `Setup ssh` once more and retry your password.

Computer
^^^^^^^^

Codes are attached to a computer, which is a set of configurations associated with a remote host machine. You can provide these configurations in **Set up a computer in AiiDA**.

.. figure:: /_static/images/setup_computer.png
   :width: 12cm

Transport & Scheduler
"""""""""""""""""""""

For `Transport` and `Scheduler`, the ``core.`` prefix is an AiiDA artifact.

The ``core.local`` transport type is set when using a local machine. As our task here pertains to a remote machine, we will be using ``core.ssh``.

Similarly, the ``core.direct`` scheduler is used for local machines. Please choose the relevant scheduler for your supercomputing center.

Prepend/Append Text
"""""""""""""""""""

AiiDA will construct the job execution command for you. However, your scheduler will need additional information regarding the job. You may add these as **Prepend text**. AiiDA will add it before the execution command.

For example, you may choose to provide the following to a `Slurm` scheduler:

.. code-block::

   #SBATCH --partition=...
   #SBATCH --constraint=...
   #SBATCH --cpus-per-task=...

   export OMP_NUM_THREADS=...
   source ...
   ulimit -s unlimited

   <execution command>

.. dropdown:: CSCS Account Specification

   CSCS machines require the user to specify the account relating to the user's project using the ``#SBATCH --account`` option. You may do so here.

You may also add commands to be run post-submission as **Append text**. AiiDA will add it after the execution command.

.. tip::

   The computer's appended text is *not* meant to be used for calculation-specific post-processing, etc., as codes come with their own appended text option (see next section). You may consider using the computer's appended text for machine-related tasks, such as clean up, file management, etc.

----

Once you've entered all the necessary details, click `Setup computer` to complete the setup. To test the new computer, click `Test computer` and wait about a minute. If all is well, you should see something similar to the following:

.. code-block::

   Report: Testing computer<daint-mc> for user<aiida@localhost>...
   * Opening connection... [OK]
   * Checking for spurious output... [OK]
   * Getting number of jobs from scheduler... [OK]: 2255 jobs found in the queue
   * Determining remote user name... [OK]: xingwang
   * Creating and deleting temporary file... [OK]
   Success: all 5 tests succeeded

Code
^^^^

Once a computer is defined, we may proceed to set up our codes in **Set up a code in AiiDA**.

.. figure:: /_static/images/setup_code.png
   :width: 12cm

For each of our codes (pw, dos, projwfc), you need to:

#. Provide a label (used in the AiiDA database)
#. Select the computer (defined in the previous section)
#. Select a code plugin (e.g. quantumespresso.pw for pw.x)
#. (**optional**) Provide a description
#. Set the absolute path to the executable on the remote machine

You may choose to prepend/append commands to the execution command that are code-specific. You may do so here. AiiDA will add it around the execution command, between the prepend/append commands that may have been added in the previous section. The final submission script would have the following form:

.. code-block::

   <computer-prepend>
   <code-prepend>
   <execution command>
   <code-append>
   <computer-append>

.. tip::

   As mentioned in the previous section, computer prepend/append text are machine-specific commands required for all jobs, while code prepend/append texts are calculation-specific. Use these texts, for example, to load necessary modules prior to execution, and to perform post-processing post-execution.

Once you are finished, click `Setup code` to complete the setup.

.. important::

   AiiDA can only make use of executables if they exist at their respective provided absolute paths. Please make sure that the codes are pre-compiled on the remote machine and are present before using them.

Once a code is created, you can close the setup widget by clicking on `Setup new code`. The newly created code is now available in its respective dropdown menu.

.. figure:: /_static/images/select_new_code.png
   :width: 12cm

----

Repeat the above process for each of the codes required for your calculation. You can skip the connection setup, as one has already been established. Similarly, you don't need to set up a new computer. But do make sure to select the same computer for each of the codes.

.. _CSCS: https://www.cscs.ch/
