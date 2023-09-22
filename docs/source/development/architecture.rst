.. _develop:architecture:

************************
Architecture
************************

Wizards UI
==========

Quantum ESPRESSO app (QeApp) uses the Wizards UI, which divides one calculation into four steps. Each step may contain several sections (panels), as shown below.

.. image:: ../_static/images/plugin_step.png

Parameter transfer between steps
---------------------------------
Data is passed to the next step by linking it to the corresponding trailet of the step. For example, the ``confirmed_structure`` of step 1 is linked to the ``input_structure`` trailet of step 2.

.. code:: python

    ipw.dlink(
            (self.structure_selection_step, "confirmed_structure"),
            (self.configure_qe_app_work_chain_step, "input_structure"),
        )

In the configuration step, there are several panels. The parameters from these panels are generated and stored as a dictionary, which is linked to the ``input_parameters`` trailet of the next submit step. The dictionary has the following structure:

.. code:: python

    {
        "workchain": {...},
        "advanced": {...},
        "bands": {...},
        "plugin_1": {...},
        ...
    }

Plugin
======
QeApp supports runing multiple properties (bands, pdos, etc) calculations in one app. Take into account the following facts:

- the configuration for a property calculation has its settings unrelated to other properties.
- the sub-workchain of the properties can be run independently.
- the analysis of the results of the properties is also independent.

Thus, we can develop a property separately and integrate it into the QeApp as a plugin. Each plugin resposible for one property calculation. For example, we could create a PDOS plugin, including its settings, workchain, and result analysis. The GUI of the PDOS plugin is only loaded when the user selects to run the PDOS property.

A QeApp plugin will typically register new panels (setting, result), and workchain to extend the app's functionality. The plugin design makes the QeApp more modularized and pluggable. So the developer can maintain their plugin as a separate folder in the QeApp (even a separate package).
