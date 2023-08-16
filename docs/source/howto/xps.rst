============================
How to run a XPS calculation
============================

Install
===========
Since the XPS module is under development. We provide a separate image for it. Set up a new `QE-XPS` profile that using the XPS image with:

   .. code-block:: console

      aiidalab-launch profile add --image superstar54/qe:xps QE-XPS

   At the prompt, enter `n` to skip editing the profile settings.


Start AiiDAlab
==============
Start AiiDAlab with the new `QE-XPS` profile:

   .. code-block:: console

      aiidalab-launch start -p QE-XPS


Run a XPS calculation for O 1s in H2O
======================================

Step 1 Select a structure
--------------------------------
For this tutorial task, please use the `From Exmple` tab, and select the `Water molecule` structure.

In the end, click the `Confirm` button.


Step 2 Configure workflow
--------------------------------

In the `Workflow` tab, set the parameters which are highlighted as follows:

- In the `Structure optimization` section, select ``Structure as is``.
- In the `properties` section, select ``XPS``


Then go to the `Basic setting` tab, and the set the parameters which are highlighted as follows:

- **Electronic Type** to `Insulator`
- In the **protocol** section, please select `fast` as the protocol for this calculation, in order to quickly get the result.




Then go to the `Advanced setting` tab, and the set the parameters which are highlighted as follows:

- In the `Pseudo Potential Selector` section, please cahnge the `Exchange-correlation functional` to `PBE`.

Then go to the `XPS setting` tab, and the set the parameters which are highlighted as follows:

- In the `Structure` section, select `molecule`.
- In the `Select element` section, select `O_1s`


In the end, click the `Confirm` button.

Step 3 Choose computational resources
---------------------------------------

We have set up a default computer (`localhost`) and the following codes:

- pw-7.2@localhost

For this calculation, you can use the default codes.

In the `Resources` section, you can set the number of CPUs to 2 for the test.

In the end, click the `Submit` button.


Step 4 Check the status and results
-----------------------------------------
The job may take ~2 minutes.


When the job is finished; then you can view the `XPS` in the `Step 4` section. If there is no result, please refresh the page, and choose the result at the top of the page in the ``Select computed workflow or start a new one``.

Here is the result:

.. figure:: /_static/images/xps_result.png
   :width: 12cm


Your own test calculation
===========================
Currently, only the XPS spectra of elements C, O, F, Si and Pt are supported.
