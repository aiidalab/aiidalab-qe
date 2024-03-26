

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
        pip: "aiidalab-qe-xyz"

3. Submit a Pull Request. Direct it to `this repository's Pull Requests section <https://github.com/aiidalab/aiidalab-qe/pulls>`_.

Plugin Entry Requirements
-------------------------

**Required Keys**

- **Top-level key:** The plugin's distribution name, which should be lowercase and prefixed by ``aiidalab-`` or ``aiida-``. For example, ``aiidalab-qe-coolfeature`` or ``aiidalab-neutron``.
- **description:** A brief description of your plugin.

**Optional Keys**

- **github:** If provided, this should be the URL to the plugin's GitHub homepage.

At least one of ``github`` or ``pip`` is required.

- **pip:** The PyPI package name for your plugin, useful for installation via pip. Example: ``aiida-quantum``.
- **documentation:** The URL to your plugin's online documentation, such as ReadTheDocs.
- **author:** The developer of the plugin.

By following these guidelines, you can ensure your plugin is correctly listed and accessible within the AiiDAlab Quantum ESPRESSO app, facilitating its discovery and use by the community.
