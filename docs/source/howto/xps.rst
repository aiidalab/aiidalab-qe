============================
How to calculate XPS spectra
============================

Overview
========
This tutorial will guide you through the process of setting up and running an XPS calculation for CH\ :sub:`3`\CHO molecule.


Steps
=====

To start, go ahead and :doc:`launch </installation/launch>` the app, then follow the steps below.


Step 1 Select a structure
--------------------------------
For this tutorial task, please use the `From Example` tab, and select the CH\ :sub:`3`\CHO molecule structure.

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
    At present, core-hole pseudopotentials for Si and O are only available for the PBE functional.

Then go to the **XPS setting** tab, and the set the parameters which are highlighted as follows:

- In the **Select core-level** section, select ``C_1s`` by ticking the appropriate box.

.. image:: ../_static/images/xps_step_2_setting_tab.png
   :align: center


Then, lick the **Confirm** button.


Step 3 Choose computational resources
---------------------------------------

For this small system, we can use the local computer (``localhost``) to run the calculation. The code ``pw-7.2@localhost`` should already be selected by default.

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

   One can read the exact binding energies from the the output of the calculation, by clicking the `outputs` tab on the node tree of the WorkChain, as shown below.

   .. figure:: /_static/images/xps_step_4_output.png
      :align: center


   When comparing the these binding energies to the experimental data, please correct the energies by the given offset listed in the ``Offset Energies`` table in the button of the XPS result tab.
   Besides, the DFT calculated binding energies do not include spin-orbit splitting of the core level state.
   We can include the spin-orbit splitting using its experimental value.
   Take `f` orbit as a example, we need subtracting :math:`3/7` of the experimental spin-orbit splitting or adding :math:`4/7` of the DFT calculated value, to get the position of the :math:`4f_{7/2}` and :math:`4f_{5/2}` peaks, respectively. Here is a table of the spin-orbit splitting for different orbitals.

   +----------------+-------------------+-------------------+
   | Orbit          | Substracting      | Adding            |
   +================+===================+===================+
   | 1s             | 0                 | 0                 |
   +----------------+-------------------+-------------------+
   | 2p             |   :math:`1/3`     |  :math:`2/3`      |
   +----------------+-------------------+-------------------+
   | 3d             | :math:`2/5`       |  :math:`3/5`      |
   +----------------+-------------------+-------------------+
   | 4f             | :math:`3/7`       |  :math:`4/7`      |
   +----------------+-------------------+-------------------+

Congratulations, you have finished this tutorial!


A real world example
====================
ETFA is commonly used as example for XPS measurements and calculations due to the extreme chemical shifts of its four different carbon atoms. [1]

.. tip::

   One can select the ETFA molecule from the `From Example` tab, and follow the same steps as above to run the XPS calculation for this molecule.
   One need use a high-performance computer to run XPS calculation for this system. Please read the relevant :doc:`How-To </howto/setup_computer_code>` section to setup code on a remote machine.


Here is the result of the XPS calculation for the ETFA molecule.


.. figure:: /_static/images/xps_etfa_dft.png
   :align: center

Here is the chemical shift from experiment. [1]

.. figure:: /_static/images/xps_etfa_exp.jpg
   :align: center


The calculated relative shifts align well with the trends observed in experimental data, underscoring the reliability of DFT calculations.
Although there are minor discrepancies in the absolute shift values, this is a recognized limitation stemming from the approximations in the exchange-correlation functional within DFT frameworks. [2]

Questions
=========

If you have any questions, please, do not hesitate to ask on the AiiDA discourse forum: https://aiida.discourse.group/.

[1] O. Travnikova *et al.*, , *Relat. Phenom.* 185, 191 (2012)
https://doi.org/10.1016/j.elspec.2012.05.009
[2] B.P. Klein  *et al.*, , *J. Phys. Condens. Matter* 33, 154005 (2021)
https://doi.org/10.1088/1361-648X/abdf00
