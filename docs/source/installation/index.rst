.. _installation:

==================
Installation Guide
==================

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   docker

***************
Access AiiDAlab
***************

As a user, you have three options to access AiiDAlab:

 1. Log into one of the `open AiiDAlab servers <https://www.aiidalab.net/deployments/>`_.
 2. Run the AiiDAlab docker container directly :ref:`on your local machine <docker>`.
 3. Download the `Quantum Mobile Virtual Machine <https://quantum-mobile.readthedocs.io/>`_, open a terminal and run ``aiidalab``.

.. _install_qeapp:

********************************
Install the Quantum-Espresso App
********************************

AiiDAlab provides for the distribution of a collection of notebooks in the form of apps through the `AiiDAlab app store`_.
These apps can be installed, updated, and removed either via the command line using the ``aiidalab`` command, or via the App store.

.. tip::

   All :fa:`terminal` **Terminal** commands shown below, can in principle also be executed directly from the host via AiiDAlab launch's ``exec`` function, for example:

   .. code-block:: console

      $ aiidalab-launch exec -- aiidalab install <app-name>

You can install new apps either using the graphical user interface via the :fa:`puzzle-piece` **App Store** or on the command line (:fa:`terminal` **Terminal**).

.. tabbed:: App Manager

    .. panels::
       :container: container-lg pb-3
       :column: col-lg-12 p-2

       **Step 1: Open the App Store.**

       Simply open AiiDAlab in the browser and click on the :fa:`puzzle-piece` icon in the top navigation bar.

       .. image:: ../_static/nav-bar-app-store.png

       This will open the app store page in a new window or tab.

       ---

       **Step 2: Search for the app you would like to install.**

       Optionally, select one or multiple categories to filter by:

       .. image:: ../_static/app-management-app-store.png

       Then scroll down until you find the app you would like to install.
       An app that is not installed yet, will be presented like this:

       .. image:: ../_static/app-management-app-not-installed.png

       Clicking on the **Install** button will install the app and its dependencies.

       In some cases the app developers will push prereleases which can be installed by clicking on the *Include prereleases* check box.
       Use this option only if you require access to a not yet released feature or you would like to test a new app version and provide feedback to the developer(s).

       ---

       **Step 3: Wait for the installation process to complete.**

       The current process for installing the app and its dependencies will be displayed via a terminal widget.
       Wait until the process has completed:

       .. image:: ../_static/app-management-app-installation-completed.png

       ---

       **Step 4: Start the app from the start page.**

       The newly installed app should now show up on the start page.

       .. image:: ../_static/app-management-start-page.png

       Each app banner also shows an indicator about whether there is an update available (see screenshot above).
       To *update the app*, click on **Manage App** and then on the **Update** buttons.


.. tabbed:: Terminal

    .. panels::
       :container: container-lg pb-3
       :column: col-lg-12 p-2

       **Step 1: Open the Terminal.**

       Open the :fa:`terminal` by clicking on the corresponding icon in the nav bar.

        .. image:: ../_static/nav-bar-terminal.png

       ---

       **Step 2: Install the app with the aiidalab command**

       .. code-block:: console

          $ aiidalab install <app-name>

       Replace ``<app-name>`` with the name of the app you would like to install, e.g., ``aiidalab install quantum-espresso``.
       Use ``aiidalab search`` to search among available apps and their versions.
       Similarly, the ``aiidalab list`` lists all currently installed apps and their versions.


.. _app-store:upgrade:

*****************************
Upgrade (or downgrade) an app
*****************************

.. tabbed:: App Manager

    .. panels::
       :container: container-lg pb-3
       :column: col-lg-12 p-2

       **Step 1: Find the app you would like to upgrade on the start page.**

       On the home app start page, simply look for the app you would like to upgrade.

       .. image:: ../_static/app-management-start-page-upgrade-available.png

       Click on the **Manage App** button to open the app manager.

       ---

       **Step 2: Open the App Management page**

       The green :fa:`arrow-circle-up` **Update** button indicates that there is a newer version of the app available.

       .. image:: ../_static/app-management-upgrade-available.png

       Click on the :fa:`arrow-circle-up` **Update** button to upgrade the app.

       By default, the app will be upgraded to the latest available version, howevever you can alternatively select any available version, including a version that is lower than the currently installed one.

.. tabbed:: Terminal

    Within the :fa:`terminal` Terminal, execute the following command to upgrade:

    .. code-block:: console

       $ aiidalab install <app-name>

    This will install the most recent version of an app, regardless of whether it is already installed or not.
    You will be prompted to confirm the operation.

    You can install a specific version, by using standard `PEP 440 version specifiers`_, for example:

    .. code-block:: console

       $ aiidalab install quantum-espresso==v22.01.0

.. _app-store:uninstall:

*****************
Uninstall an app
*****************

Uninstalling an app works similar to upgrading or downgrading an app via the **Manage App** page or on the Terminal.

.. tabbed:: App Manager

    .. panels::
       :container: container-lg pb-3
       :column: col-lg-12 p-2

       **Step 1: Find the app you would like to uninstall on the start page.**

       On the home app start page, simply look for the app you would like to uninstall.

       .. image:: ../_static/app-management-start-page.png

       Click on the **Manage App** button to open the app manager.

       ---

       **Step 2: Uninstall**

       The app manager allows you to uninstall the app or to install a different version.

       .. image:: ../_static/app-management-app-installed.png

       Click on the :fa:`trash` **Uninstall** button to uninstall the app.

       .. note::

          In some cases you will see a warning that uninstalling the app might lead to data loss.
          That warning indicates that there are local modifications to the app source code.
          You can safely ignore this warning and click on the "Ignore" check box in case that you are sure that any local modifications are safe to delete.

.. tabbed:: Terminal

    Within the :fa:`terminal` Terminal, execute the following command to uninstall an app:

    .. code-block:: console

       $ aiidalab uninstall <app-name>

    You will be prompted to confirm the operation.

.. _AiiDAlab app store: https://aiidalab.github.io/aiidalab-registry
.. _PEP 440 version specifiers: https://www.python.org/dev/peps/pep-0440/#version-specifiers
