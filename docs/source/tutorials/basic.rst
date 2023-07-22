.. _run_calculation_on_local_computer:

Run calculation on local computer
=================================

Goal: running a band structure calculation on Silicon using the AiiDAlab machine. Get the band structure and PDOS.


As shown at the top of the page, there will be four steps to carry out a QE calculation.

Step 1 Select a structure
-------------------------

You can select a structure from four sources.

- Upload file
- OPTIMADE
- AiiDA database
- From Examples

For the first task, please use the `From Examples` tab, and select `Silicon (diamond)`.

.. figure:: /_static/images/step1_select_structure.png
   :width: 12cm


In the end, click the `Confirm` button.



Step 2 Configure workflow
-------------------------

Set the input parameters

- There are three relax types:

    - Structure as is
    - Atomic positions
    - Full geometry

- There are two properties that can be computed:
    - Band structure
    - Projected density of states

For this test calculation, please click all of them.

Please select `fast` as the protocol for this calculation, in order to quickly get the result, since you only have two CPUs both for running AiiDAlab and to run the QE calculation.
Note, of course, that this is only for demo purposes: for any production calculation, please set up a "real" remote computer (cluster, supercomputer) and use it!
You can check the next section on how to setup a remote computer, and the Quantum ESPRESSO codes to run a band structure calculation on it.

.. figure:: /_static/images/step2_configure.png
   :width: 12cm

In the end, click the `Confirm` button.


Step 3 Choose computational resources.
--------------------------------------

**Codes**

We have set up a default computer (`localhost`) and the following codes:

- pw-7.0@localhost
- dos-7.0@localhost
- projwfc-7.0@localhost

For this calculation, you can use the default codes.

**Resources**
You can set the number of CPUs to 1 for the test.


.. figure:: /_static/images/step3_resources.png
   :width: 12cm


Click the `Submit` button and the calcuation is launched.

.. note::

    There is no need to wait for the computation to finish: you can go head and submit a new calculation in parallel.

Step 4 Check the status and results.
------------------------------------

The job may take ~3 minutes. While waiting, you can also check the job information by clicking the Workchain item in the tree view.

In addition, above `Step 1`, there is a textbox showing the status of a workflow.


.. figure:: /_static/images/workchain_selector.png
   :width: 12cm


When the job is finished, refresh the if needed; then you can view the `Final Geometry` and the `Electronic Structure` in the `Step 4` section.

.. figure:: /_static/images/step4_results.png
   :width: 12cm
