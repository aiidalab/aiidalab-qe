.. _install:

Install the QE app
==================

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

       .. image:: ../_static/images/nav-bar-app-store.png

       This will open the app store page in a new window or tab.

       ---

       **Step 2: Search for the app you would like to install.**

       Optionally, select one or multiple categories to filter by:

       .. image:: ../_static/images/app-management-app-store.png

       Then scroll down until you find the app you would like to install.
       An app that is not installed yet, will be presented like this:

       .. image:: ../_static/images/app-management-app-not-installed.png

       Clicking on the **Install** button will install the app and its dependencies.

       In some cases the app developers will push prereleases which can be installed by clicking on the *Include prereleases* check box.
       Use this option only if you require access to a not yet released feature or you would like to test a new app version and provide feedback to the developer(s).

       ---

       **Step 3: Wait for the installation process to complete.**

       The current process for installing the app and its dependencies will be displayed via a terminal widget.
       Wait until the process has completed:

       .. image:: ../_static/images/app-management-app-installation-completed.png

       ---

       **Step 4: Start the app from the start page.**

       The newly installed app should now show up on the start page.

       .. image:: ../_static/images/app-management-start-page.png

       Each app banner also shows an indicator about whether there is an update available (see screenshot above).
       To *update the app*, click on **Manage App** and then on the **Update** buttons.


.. tabbed:: Terminal

    .. panels::
       :container: container-lg pb-3
       :column: col-lg-12 p-2

       **Step 1: Open the Terminal.**

       Open the :fa:`terminal` by clicking on the corresponding icon in the nav bar.

        .. image:: ../_static/images/nav-bar-terminal.png

       ---

       **Step 2: Install the app with the aiidalab command**

       .. code-block:: console

          $ aiidalab install <app-name>

       Replace ``<app-name>`` with the name of the app you would like to install, e.g., ``aiidalab install quantum-espresso``.
       Use ``aiidalab search`` to search among available apps and their versions.
       Similarly, the ``aiidalab list`` lists all currently installed apps and their versions.

.. _AiiDAlab app store: https://aiidalab.github.io/aiidalab-registry
