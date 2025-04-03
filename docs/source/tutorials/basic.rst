==============
Basic Tutorial
==============

In this tutorial, we will show how to use the basic features of the Quantum ESPRESSO app to run a geometry optimization of bulk Silicon and obtain its band structure.

----

The Quantum ESPRESSO app involves the following steps:

#. :ref:`Selecting a structure <basic:structure>`
#. :ref:`Configuring the workflow <basic:config>`
#. :ref:`Choosing computational resources and submitting the job <basic:resources>`
#. :ref:`Monitoring the status and retrieving results <basic:results>`

The sections below provide a brief overview of each step, along with screenshots from the app.

----

.. _basic:structure:

Structure selection
*******************

You can select a structure from the AiiDA database, upload a custom structure, use the OPTIMADE service to search for a structure across registered databases, or choose from a list of examples.
Once uploaded, you can visualize and modify the structure using the built-in visualization tool.

.. figure:: ../_static/images/in_app_guides/structure_selection.png
   :width: 100%
   :align: center
   :class: img-responsive

.. _basic:config:

Workflow configuration
**********************

In the configuration step, you define the workflow by including relaxation (optional) and selecting the desired properties to compute (e.g., band structure, density of states, etc.), as well as specifying parameters for the selected calculations (e.g., protocol, k-point grid, cutoffs, etc.).

.. figure:: ../_static/images/in_app_guides/workflow_configuration.png
   :width: 100%
   :align: center
   :class: img-responsive

.. _basic:resources:

Resources and submission
************************

In this step, you select the computer and code to use for the calculation, as well as specify the number of nodes and CPUs. When you are ready, you can submit the job.

.. figure:: ../_static/images/in_app_guides/computational_resources.png
   :width: 100%
   :align: center
   :class: img-responsive

.. _basic:results:

Monitoring and results
**********************

Once submitted, the app redirects you to the last step, from where you can view a summary of parameters, monitor the status of the workflow, and visualize the results of each calculation as they become available.

.. figure:: ../_static/images/in_app_guides/results.png
   :width: 100%
   :align: center
   :class: img-responsive
