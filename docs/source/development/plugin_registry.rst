

Plugin Registry
=========================================

If you are either in the process of creating a new plugin or already have one developed, you're encouraged to register your plugin here to become part of the official AiiDAlab Quantum ESPRESSO App plugin ecosystem.

Registering Your Plugin
-----------------------

To include your plugin in the registry, follow these steps:

1. Fork this `repository <https://github.com/aiidalab/aiidalab-qe>`_.

2. Add your plugin to the `plugins.yaml` file. Place your entry at the end of the file, following this example:

   .. code-block:: yaml

      aiidalab-qe-xyz:
        description: "Quantum ESPRESSO plugin for XYZ by AiiDAlab."
        author: "Alice Doe"
        github: "https://github.com/alicedoe/aiidalab-qe-xyz"
        documentation: "https://aiidalab-qe-xyz.readthedocs.io/"
        pip: "aiidalab-qe-xyz==version-of-the-code"
        post-install: "post-install-command"

3. Submit a Pull Request. Direct it to `this repository's Pull Requests section <https://github.com/aiidalab/aiidalab-qe/pulls>`_.

Plugin Entry Requirements
-------------------------

**Required Keys**

- **Top-level key:** The plugin's distribution name, which should be lowercase and prefixed by ``aiidalab-`` or ``aiida-``. For example, ``aiidalab-qe-coolfeature`` or ``aiidalab-neutron``.
- **description:** A brief description of your plugin.

**Optional Keys**

- **github:** If provided, this should be the URL to the plugin's GitHub homepage.

At least one of ``github`` or ``pip`` is required. ``pip`` installation will be preferred if both are provided, and "==version-of-the-code" can be omitted (but strongly suggested, to ensure compatiblity).

- **pip:** The PyPI package name for your plugin, useful for installation via pip. Example: ``aiida-quantum``.
- **documentation:** The URL to your plugin's online documentation, such as ReadTheDocs.
- **author:** The developer of the plugin.
- **post-install:** a post install Command Line Interface (CLI) command which should be defined inside your plugin if you needs it. For example in the ``aiidalab-qe-vibroscopy`` plugin, we automatically setup the phonopy code via this command. See below for more explanations.

How to define a post install command in your plugin
---------------------------------------------------------------------

If you need to run a post-install command, you can define it in the CLI of your package, meant to be run as ``package-name post-install-command``.
You can define the CLI  in the ``__main__.py``  (in your source folder) and ``pyproject.toml`` files. Have a look at the  ``aiidalab-qe-vibroscopy`` plugin for an example on how to do this.
There, the automatic setup the phonopy code in the AiiDA database is implemented.
The command will be then triggered after the installation of the plugin (only) from the plugin list page of the Quantum ESPRESSO app.

By following these guidelines, you can ensure your plugin is correctly listed and accessible within the AiiDAlab Quantum ESPRESSO app, facilitating its discovery and use by the community.
