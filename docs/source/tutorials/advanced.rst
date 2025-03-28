=================
Advanced Tutorial
=================

Overview
--------

Calculations probing more intricate properties often require additional computational resources. In this tutorial, we will again compute a band structure, but this time, we will do it for a magnetic system, namely bulk nickel. As such, will set up our environment to submit to a remote machine capable of handling the calculation.

Start
-----

.. admonition:: Goal

   To submit a band structure calculation for magnetic nickel on a remote machine.

To start, go ahead and :doc:`launch </installation/launch>` the app, then follow the steps below.

.. note::

   For brevity, we did not include images for steps 1-3. Please refer to the :doc:`basic tutorial </tutorials/basic>` for visual aid.

Step 1: Select a structure
**************************

Select `Nickel (hcp)` from the `From examples` tab and click `Confirm`.

Step 2: Configure the workflow
******************************

Select `Full geometry` to relax the structure. As we are looking at magnetic nickel, we also want to turn on `Magnetism`. Go ahead and do so now.

Select both `Band structure` and `Projected density of states` as the properties of interest.

For the protocol, select `fast` to quickly produce results, or, if you have enough resources, you can select the `stringent` protocol for increased accuracy.

Click `Confirm` to proceed.

Step 3 - Choose computational resources
***************************************

As mentioned above, our calculation will require more computational resources than the basic tutorial. As such, we will take a brief detour through the relevant :doc:`How-To </howto/setup_computer_code>` section to configure the environment with a remote machine. When you're done, click `Submit` to submit the calculation.

Step 4: Check the status and results
************************************

Again, while the calculation is running, you can monitor its status as shown in the :ref:`basic tutorial <basic_status>`. You can view the results once the calculation is finished.

Summary
-------

In this tutorial, you learned how to submit a band structure calculation on a remote machine using the Quantum ESPRESSO app.
