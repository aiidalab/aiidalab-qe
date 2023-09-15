.. _develop:create-plugin:

************************
Create your plugin
************************

Quantum ESPRESSO app (QeApp) uses the Wizards UI, which divides one calculation into four steps. Each step may contain several sections (panels), as shown below.

.. image:: ../_static/images/plugin_step.png

A QeApp plugin will typically register new panels (setting, result), and workchain to extend the app's functionality. The plugin design make the QeApp more modularized, and pluggable. So the developer can maintain their plugin as a separate folder in the QeApp (even a separate pakcage).

Your First Add-on
================================

Here is the simplest plugin to print the formula of the input structure:

**Outline**, it will be shown as a checkbox in the workflow panel, as shown below.

.. code-block:: python

    from aiidalab_qe.common.panel import OutlinePanel

    class Outline(OutlinePanel):
        title = "Hello World"


.. image:: ../_static/images/plugin_outline.png




**Setting**, it will register a new panel in the configuration Step. In this class, one needs to implement the `get_panel_value` and `load_panel_value` methods to tell QeApp how to handle the panel values.

.. image:: ../_static/images/plugin_setting.png


.. code-block:: python

    from aiidalab_qe.common.panel import Panel

    class Setting(Panel):
        """"""
        title = "Hello world"

        def __init__(self, **kwargs):
            self.name = ipw.Text(value="", description="Your name:")
            self.children = [self.name]
            super().__init__(**kwargs)

        def get_panel_value(self):
            return {"name": Str(self.name.value)}

        def load_panel_value(self, input_dict):
            self.name.value = input_dict.get("name", 1)

**Result**, it will register a new panel in the Results Step. In this class, one needs to implement the `_update_view` method to tell QeApp how to show the outputs of the workchain.

.. image:: ../_static/images/plugin_result.png


.. code-block:: python

    from aiidalab_qe.common.panel import ResultPanel

    class Result(ResultPanel):
        title = "Hello world"
        workchain_label = "hello_world"

        def _update_view(self):
            name = self.node.outputs.name.value
            formula = self.node.outputs.structure.get_formula()
            self.summary_view = ipw.HTML(
                f"""<div> <h4>Hello {name}</h4> The input structure is: {formula} </div>""".format()
            )
            self.children = [ipw.HBox(children=[self.summary_view])]


**WorkChain**, it will be loaded into the QeApp workchain. One needs to implement a `get_builder` function to extract the input parameters for the sub-workchain.

.. code-block:: python

    def get_builder(codes, structure, parameters):
        """Get the workchain specific parameters
        """
        parameters = parameters.get("hello_world", {})
        builder = HelloWorldWorkChain.get_builder_from_protocol(
                    codes=codes,
                    structure=structure,
                    parameters=parameters,
                )
        return builder


    workchain_and_builder = {
        "workchain": HelloWorldWorkChain,
        "get_builder": get_builder,
        }

**Entry point**, here is the entry point for this plugin. One needs to add it to `entry_points` inside the setup file.

.. code-block:: python

    hello_world ={
    "outline": Outline,
    "setting": Setting,
    "workchain": workchain_and_builder,
    "result": Result,
    }


.. code-block:: python

    entry_points={
            "aiidalab_qe.properties": [
                "hello_world = aiidalab_qe_hello_world:hello_world",
            ],
        },

Note: one plugin does not need to register all the items (settings, workchain, results). The panel in each step is pluggable, which means you could only register one item in a plugin. For example, you can only add a new `Structure` panel in Step 1 without doing any property calculation.

You can add this plugin as a folder in the QeApp package, or create a new package for it.

**Bringing It All Together**, You can find all the code above in this [github repository](https://github.com/superstar54/aiidalab-qe-hello-world).
