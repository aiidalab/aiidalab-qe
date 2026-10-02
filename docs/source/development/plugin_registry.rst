

Plugin Registry
=========================================

If you are either in the process of creating a new plugin or already have one developed, you're encouraged to register your plugin here to become part of the official AiiDAlab Quantum ESPRESSO App plugin ecosystem.

Registering Your Plugin
-----------------------

To include your plugin in the registry, follow these steps:

1. Fork this `repository <https://github.com/aiidalab/aiidalab-qe>`_.

2. Add your plugin to the `plugins.yaml` file. Place your entry at the end of the file, following this example:

   .. code-block:: yaml

      Top-level key:
        title: "Description to show on top"
        description: "Quantum ESPRESSO plugin for XYZ by AiiDAlab."
        author: "Alice Doe"
        github: "https://github.com/alicedoe/aiidalab-qe-xyz"
        documentation: "https://aiidalab-qe-xyz.readthedocs.io/"
        pip: "aiidalab-qe-xyz==version-of-the-code"
        post_install: "post-install-command"

3. Submit a Pull Request. Direct it to `this repository's Pull Requests section <https://github.com/aiidalab/aiidalab-qe/pulls>`_.

Plugin Entry Requirements
-------------------------

**Required Keys**

- **Top-level key:**   The plugin's distribution name, which should be lowercase and prefixed by ``aiidalab-`` or ``aiida-``. For example, ``aiidalab-qe-coolfeature`` or ``aiidalab-neutron``.
- **title:** Brief title to show on top of the plugin entry. Should contain the main properties we can compute with the given plugin.
- **description:** A brief description of your plugin. Can include more verbose informations with respect to the title.

**Optional Keys**

- **github:** If provided, this should be the URL to the plugin's GitHub homepage.

At least one of ``github`` or ``pip`` is required. ``pip`` installation will be preferred if both are provided. Its value can include a PEP 508 version constraint, for example ``aiidalab-qe-xyz>=1.2.0`` to specify a minimum version or ``aiidalab-qe-xyz==1.2.0`` to pin an exact version. Installed plugins that do not meet their constraint are shown as incompatible in Step 2.1 and can be updated in the Plugin Store.

- **pip:** The PyPI package name for your plugin, useful for installation via pip. Example: ``aiida-quantum``.
- **documentation:** The URL to your plugin's online documentation, such as ReadTheDocs.
- **author:** The developer of the plugin.
- **post_install:** a post install Command Line Interface (CLI) command which should be defined inside your plugin if you needs it. For example in the ``aiidalab-qe-vibroscopy`` plugin, we automatically setup the phonopy code via this command. See below for more explanations.
- **requires_aiidalab_qe:** an optional version specifier for the minimum compatible AiiDAlab QE app version, such as ``>26.09.0``.

Registry Validation
-------------------

The registry is validated when it is loaded. Each entry must provide non-empty
``title`` and ``description`` values, and at least one non-empty ``pip`` or
``github`` installation source. A ``pip`` value must be a valid PEP 508
requirement, and ``requires_aiidalab_qe`` must be a valid version specifier.
Unknown keys are allowed so registry metadata can be extended without changing
the loader. The test suite also validates every entry in the checked-in
``plugins.yaml`` file.

How to define a post install command in your plugin
---------------------------------------------------------------------
If you need to run a post-install command, you can define it in the CLI of your package. The command should be designed to be run as ``package-name post-install-command``.
To define the CLI, you can use the ``__main__.py`` file in your source folder and the ``pyproject.toml`` file. You can refer to the `aiidalab-qe-vibroscopy <https://github.com/mikibonacci/aiidalab-qe-vibroscopy>`_ plugin for an example of how to do this.
In that plugin, the automatic setup for the phonopy code is implemented. It assumes that the ``phonopy`` binary is already present on the machine, as the plugin will install it as a dependency.
The post-install command is triggered after installing the plugin from the Plugin Store. If needed, it can also be run on its own from the plugin entry in the store, without reinstalling the package.
