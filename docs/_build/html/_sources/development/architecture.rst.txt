.. _develop:architecture:

************************
Architecture
************************

Wizards-like UI
================

QuantumESPRESSO app uses the `Wizards-like UI design <https://en.wikipedia.org/wiki/Wizard_(software)>`_, which divides one calculation into four steps.
Each step may contain several sections (panels), as shown in the figure below.

.. image:: ../_static/images/plugin_step.png

Parameter transfer between steps
---------------------------------

The data is passed to the next step by linking it to the corresponding ``trait`` of the step widget.
For example, the ``confirmed_structure`` of step 1 links to the ``input_structure`` trait of step 2.

.. code:: python

    ipw.dlink(
            (self.structure_selection_step, "confirmed_structure"),
            (self.configure_qe_app_work_chain_step, "input_structure"),
        )

Each step widget can contains several panels. In the configuration step, the parameters from the panels are generated and stored as a dictionary, which is linked to the ``input_parameters`` trait of the next submit step.
The dictionary has the following structure:

.. code:: python

    {
        "workchain": {
            "protocol": "fast",
            "relax_type": "positions",
            "properties": ["bands", "pdos", "relax"],
            "spin_type": "none",
            "electronic_type": "insulator",
        },
        "advanced": {
            "initial_magnetic_moments": None,
            "pw": {
                "parameters": {
                    "SYSTEM": {"ecutwfc": 30.0, "ecutrho": 240.0, "tot_charge": 0.0}
                },
                "pseudos": {"Si": "eaef3352-2b0e-4205-b404-e6565a88aec8"},
            },
            "pseudo_family": "SSSP/1.3/PBEsol/efficiency",
            "kpoints_distance": 0.5,
        },
        "bands": {},
        "pdos": {...},
        "plugin_1": {...},
        "plugin_2": {...},
    }

Plugin
=========
QuantumESPRESSO app supports running multiple properties (bands, pdos, etc.) calculations in one app.
The individual properties to be developed and seamlessly integrated into the QuantumESPRESSO app as plugins. This integration is made possible due to several key aspects:

- the configuration for a property calculation has its settings unrelated to other properties.
- the sub-workchain of the properties can be run independently.
- the analysis of the results of the properties is independent.

Each plugin responsible for one property calculation.
For instance, we could create a PDOS plugin, including its settings, workchain, and result analysis.
The GUI of the PDOS plugin is only loaded when the user selects to run it.
Here is an example, where two new setting panels are shown when the user selects to run the properties.

.. figure:: ../_static/images/plugin_example.gif

A QuantumESPRESSO app plugin will typically register new panels (setting, result), and workchain to extend the functionality of the app.
The plugin design makes the QuantumESPRESSO app more modularized and pluggable.
Consequently, developers have the flexibility to manage their plugins in a distinct folder within the QuantumESPRESSO application's codebase, or they may choose to maintain it as an independent package.
