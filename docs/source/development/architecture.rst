.. _develop:create-plugin:

************************
Architecture
************************

Wizards UI
==========

Quantum ESPRESSO app (QeApp) uses the Wizards UI, which divides one calculation into four steps. Each step may contain several sections (panels), as shown below.

.. image:: ../_static/images/plugin_step.png

Parameter trasfer between steps
-------------------------------
In each step, user can set the parameters in the panels. The parameters are generated and stored as a value (e.g. Structure or dictionary), which is passed to the next step by linking it to the cooresponding trailet of the step. For example, the ``confirmed_structure`` of step 1 is linked to the ``input_structure`` trailet of step 2.

.. code:: python

    ipw.dlink(
            (self.structure_selection_step, "confirmed_structure"),
            (self.submit_qe_app_work_chain_step, "input_structure"),
        )


Plugin
======
QeApp supports run multiple properties (bands, pdos, etc) calculations in one app. Take into account the following facts:

- the configuration for a property calculation has its settings unrelated to other properties.
- the sub-workchain of the properties can be run independently.
- the analysis of the results of the properties is also independent.

Thus, we can divide the App into several small groups, and each group is responsible for one property calculation. For example, we could create a group for the `PDOS` property, which includes its settings, workchain, and result analysis. The `PDOS` group is only loaded when the user selects to run the `PDOS` property

A QeApp plugin will typically register new panels (setting, result), and workchain to extend the app's functionality. The plugin design make the QeApp more modularized, and pluggable. So the developer can maintain their plugin as a separate folder in the QeApp (even a separate pakcage).
