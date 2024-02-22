============================
How to calculate XPS spectra
============================

Overview
========
This tutorial will guide you through the process of setting up and running a XPS calculation for CH\ :sub:`3`\CHO molecule.


Steps
=====

To start, go ahead and :doc:`launch </installation/launch>` the app, then follow the steps below.


Step 1 Select a structure
--------------------------------
For this tutorial task, please use the `From Exmple` tab, and select the CH\ :sub:`3`\CHO molecule structure.

In the end, click the `Confirm` button.

.. figure:: /_static/images/xps_step_1.png
   :align: center


Step 2 Configure workflow
--------------------------------

In the **Basic Setting** tab, set the parameters which are highlighted as follows:

- In the **Structure optimization** section, select ``Structure as is``.
- Set **Electronic Type** to ``Insulator``
- In the **properties** section, select ``X-ray photoelectron spectroscopy (XPS)``
- In the **protocol** section, please select `fast` as the protocol for this calculation, in order to quickly get the result.


Then go to the **Advanced settings** tab, navigate to `Accuracy and precision`, tick the `Override` box on the right-hand-side and in the dropdown box under `Exchange-correlation functional` select `PBE`.

.. image:: ../_static/images/XAS_Plugin-Set_Adv_Options-Alt-Annotated-Cropped.png
   :align: center


.. note::
    At present, the core-hole pseudopotential for Si and O only available in PBE functional.

Then go to the **XPS setting** tab, and the set the parameters which are highlighted as follows:

- In the **Select core-level** section, select ``C_1s``.


Then, lick the **Confirm** button.


Step 3 Choose computational resources
---------------------------------------

For this small system, we can use the local computer (``localhost``) to run the calculation. The code ``pw-7.2@localhost`` is already selected by default.

In the **Resources** section, you can set the number of CPUs to 4 for the test.

Then, click the **Submit** button.

.. note::

   In pratice, one usually need high-performance computer to run XPS calculation for large system. Please read the relevant :doc:`How-To </howto/setup_computer_code>` section to setup code on a remote machine.




Step 4 Check the status and results
-----------------------------------------
The job may take ~2 minutes to finish.
While the calculation is running, you can monitor its status as shown in the :ref:`basic tutorial <basic_status>`.
When the job is finished, you can view result spectra in the `XPS` tab.

.. tip::

   If the `XPS` tab is now shown when the jobs is finished.
   Click the ``QeAppWorkChain<pk>`` item on top of the nodetree to refresh the step.

Here is the result of the XPS calculation.
You can click the **Binding energy** button to view the calculated binding energies.
You can change which element to view XPS spectra for using the dropdown box in the top left.

.. figure:: /_static/images/xps_step_4.png
   :align: center

.. tip::

   When comparing the calculated binding energies to the experimental data, please correct the energies by the given offset listed in the ``Offset Energies`` table.
   Besides, the DFT calculated binding energies do not include spin-orbit splitting of the core level state.
   We can include the spin-orbit splitting using its experimental value.
   Take `f` orbit as a example, we need subtracting :math:`3/7` of the experimental spin-orbit splitting or adding :math:`4/7` of the DFT calculated value, to get the position of the :math:`4f_{7/2}` and :math:`4f_{5/2}` peaks, respectively.


Congratulations, you have finished this tutorial!


A real world example
====================

Here is a XPS calculation for ETFA molecule, which is also referred to as the “ESCA molecule” in the literature. [1]

.. figure:: /_static/images/xps_etfa_dft.png
   :align: center

One can compared the calculated chemical shift with the experimental data.

.. figure:: /_static/images/xps_etfa_exp.jpg
   :align: center


Questions
=========

If you have any questions, please, do not hesitate to ask on the AiiDA discourse forum: https://aiida.discourse.group/.

[1] Travnikova, O. et al., Relat. Phenom. 2012, 185, 191-197
https://doi.org/10.1016/j.elspec.2012.05.009
