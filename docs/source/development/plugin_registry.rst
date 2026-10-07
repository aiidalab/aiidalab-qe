

Plugin Registry
=========================================

If you are developing an AiiDAlab Quantum ESPRESSO plugin, you can register it in the app's plugin catalog.

Registering Your Plugin
-----------------------

To include your plugin in the catalog, follow these steps:

1. Fork this `repository <https://github.com/aiidalab/aiidalab-qe>`_.

2. Add your plugin to ``src/aiidalab_qe/plugins/plugins.yaml``. Place your entry at the end of the file, following this example:

   .. code-block:: yaml

      Top-level key:
         title: "XYZ"
         description: "Quantum ESPRESSO plugin for XYZ."
         pip: "aiidalab-qe-xyz>=1.2,<2"
         github: "https://github.com/alicedoe/aiidalab-qe-xyz"
         author: "Alice Doe"
         documentation: "https://aiidalab-qe-xyz.readthedocs.io/"
         post_install: "post-install-command"
         status: "stable"
         category: "calculation"
         requires_aiidalab_qe: ">=26.09,<27"

3. Submit a Pull Request. Direct it to `this repository's Pull Requests section <https://github.com/aiidalab/aiidalab-qe/pulls>`_.

Plugin Entry Requirements
-------------------------

The catalog is bundled with the app release at ``src/aiidalab_qe/plugins/plugins.yaml``. Compatibility requirements therefore describe the plugin versions supported by that app release; update the catalog through a pull request when support changes.

**Required Keys**

- **Top-level key:** A unique registry identifier and installed distribution name, conventionally lowercase (for example, ``aiidalab-qe-xyz``).
- **title:** Brief title to show on top of the plugin entry. Should contain the main properties we can compute with the given plugin.
- **description:** A brief description of your plugin.

**Optional Keys**

- **pip:** PEP 508 package requirement used to install the plugin. If both ``pip`` and ``github`` are provided, ``pip`` is used.
- **github:** GitHub repository URL used as an install source when ``pip`` is omitted. A ref may be specified with ``@ref``.
- **author:** Plugin authors.
- **documentation:** URL to the plugin documentation.
- **post_install:** Optional CLI command to run after installation or update. The plugin must expose it as ``python -m <package_module> <command>``; it should be safe to rerun.
- **status:** Display label such as ``experimental``, ``beta``, ``stable``, ``production``, ``deprecated`` or ``archived``.
- **category:** Plugin category. Use ``calculation`` for a property selectable in Step 2.1; other categories are managed by the Plugin Manager but are not listed as calculation properties.
- **requires_aiidalab_qe:** Optional PEP 440 version specifier for compatible app versions.

At least one of ``pip`` or ``github`` is required. Use a PEP 508 requirement in ``pip`` to constrain the plugin version. For example, ``aiidalab-qe-xyz>=1.2,<2`` accepts compatible 1.x releases, while ``aiidalab-qe-xyz==1.2.3`` pins one exact version. Quote the whole requirement in YAML when it contains special characters or commas.

Use ``requires_aiidalab_qe`` to constrain the app version independently of the plugin package version. It accepts PEP 440 specifiers, including ranges such as ``>=26.09,<27`` or ``>=26.09.1,!=26.10.0``. A value like ``>=26.09`` is a minimum, not an upper bound; include an upper bound when compatibility is limited to a release series.

The Plugin Manager evaluates the installed distribution version against the ``pip`` requirement and the installed app version against ``requires_aiidalab_qe``. In a new calculation, incompatible plugins are excluded from Step 2.1 and cannot be selected. When an existing process is loaded, its saved selections remain visible in Steps 2-3; Step 4 independently blocks incompatible result panels and recommends opening the Plugin Manager. Update runs the configured post-install command and plugin-loading validation before restarting the daemon.

Registry Validation
-------------------

The registry is validated when it is loaded. Each entry must provide non-empty
``title`` and ``description`` values, and at least one non-empty ``pip`` or
``github`` installation source. A ``pip`` value must be a valid PEP 508
requirement, and ``requires_aiidalab_qe`` must be a valid version specifier.
Unknown keys are allowed so registry metadata can be extended without changing
the loader. The test suite also validates every entry in the checked-in
``src/aiidalab_qe/plugins/plugins.yaml`` file.

How to define a post install command in your plugin
---------------------------------------------------------------------
If you need to run a post-install command, you can define it in the CLI of your plugin. The command should be designed to be run as ``plugin-name post-install-command``.
To define the CLI, you can use the ``__main__.py`` file in your source folder and the ``pyproject.toml`` file. You can refer to the `aiidalab-qe-vibroscopy <https://github.com/mikibonacci/aiidalab-qe-vibroscopy>`_ plugin for an example of how to do this.
In that plugin, the automatic setup for the phonopy code is implemented. It assumes that the ``phonopy`` binary is already present on the machine, as the plugin will install it as a dependency.
The post-install command is triggered after installing the plugin from the Plugin Store. If needed, it can also be run on its own from the plugin entry in the store, without reinstalling the plugin.
