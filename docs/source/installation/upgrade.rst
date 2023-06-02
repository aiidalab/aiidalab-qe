.. _upgrade:

===============
Upgrade the app
===============

.. tabbed:: App Manager

    .. panels::
       :container: container-lg pb-3
       :column: col-lg-12 p-2

       **Step 1: Find the app you would like to upgrade on the start page.**

       On the home app start page, simply look for the app you would like to upgrade.

       .. image:: ../_static/images/app-management-start-page-upgrade-available.png

       Click on the **Manage App** button to open the app manager.

       ---

       **Step 2: Open the App Management page**

       The green :fa:`arrow-circle-up` **Update** button indicates that there is a newer version of the app available.

       .. image:: ../_static/images/app-management-upgrade-available.png

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

.. _PEP 440 version specifiers: https://www.python.org/dev/peps/pep-0440/#version-specifiers
