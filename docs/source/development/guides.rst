.. _develop:guides:

******
Guides
******

Throughout the app, instances of ``aiidalab_qe.common.infobox.InAppGuide`` have been dropped to provide users with helpful information in the form of toggleable tutorials.
These guide section widgets observe a global instance of ``aiidalab_qe.common.guide_manager.GuideManager`` to obtain their content and toggle their visibility.
The guide manager is initialized at app startup and scans the local codebase, as well as entry-point-based plugins, for defined guides.
The discovered guides are then made available as a nested radio button widget of the form ``category`` -> ``guide``.
This allows both built-in ``general/<guide-name>`` guides and external ``<plugin-identifier>/<guide-name>`` guides, ensuring guide uniqueness.

Developers have two options when using the ``InAppGuide`` widget:

#. **Widget-based guide sections**: Developers can define the guide content directly in the widget's constructor using ``ipywidgets`` widgets.
   When using this approach, the developer may provide a ``guide_id`` in the constructor to associate the guide section with a specific guide.
   Otherwise, if no ``guide_id`` is provided, the guide section will be displayed for any active guide.

   .. code:: python

       InAppGuide(
           guide_id="general/basic",
           children=[
               ipw.HTML("Some basic information."),
           ],
       )

#. **HTML-based guide sections**: Developers can define the guide content in an HTML document, tagging the various sections using the HTML ``id`` attribute.
   The guide manager will load the HTML document associated with the selected guide and parse it using ``BeautifulSoup``.

   .. code:: html

       <div id="structure-step">
           <p>Some basic information regarding the structure.</p>
       </div>

       <div id="configuration-step">
           <p>Some basic information regarding the configuration.</p>
       </div>

       ...

   HTML-based ``InAppGuide`` widgets are instantiated without children (or a guide id), but instead with an ``identifier`` corresponding to the relevant HTML ``id`` attribute.

   .. code:: python

       InAppGuide(identifier="structure-step")

A decent amount of ``InAppGuide(identifier=<identifier>)`` instances have been placed strategically throughout the app.
Developers may suggest additional core guide sections via GitHub pull requests.
For plugin developers, additional instances of either flavor are recommended to be added in any component of the plugin in conjunction with dedicated plugin-specific guides.

Plugin guides
-------------

Plugin developers can enhance user experience while using the app by introducing custom guides.
To do so, add the following key/value entry in your plugin's ``__init__.py`` file:

.. code:: python

    my_plugin = {
        ...
        "guides": {
            "title": <title>,
            "path": <path-to-guide>,
        },
    }

where ``title`` is the title of the guide displayed in the app, and ``path-to-guide`` is the path (``Path`` object or absolute string path) to the directory containing the guide HTML files.
On app start, the guide manager will scan the plugin entry points for the ``guides`` key and load the guides accordingly.

It is possible also to provide just the path to the guide HTML files, in which case the title will be set to the plugin identifier. In this case, just provide:

.. code:: python

    my_plugin = {
        ...
        "guides": <path-to-guide>,
    }

Guide order
-----------

When naming your guide HTML documents, prefix the file name with ``#_``. The number ``#`` will determine the order in which the guides are displayed in the list.
