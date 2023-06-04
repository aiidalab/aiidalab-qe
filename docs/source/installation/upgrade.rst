.. _upgrade:

Upgrade the app
===============

.. tabbed:: App Manager

    .. panels::
       :container: container-lg pb-3
       :column: col-lg-12 p-2

       **Step 1: Open the app management page**

       Click on the **Manage App** button to open the app manager.

       .. image:: ../_static/images/app-management-start-page-upgrade-available.png

       ---

       **Step 2: Upgrade the app**

       The green :fa:`arrow-circle-up` **Update** button indicates that a newer version is available. Click it to upgrade the app.

       .. image:: ../_static/images/app-management-upgrade-available.png

       .. note::

          By default, the app will be upgraded to the latest available version. Alternatively, you can select any available version including versions lower than the currently installed one.

.. tabbed:: Terminal

    .. panels::
       :container: container-lg pb-3
       :column: col-lg-12 p-2

       You can upgrade the app from the built-in terminal with

       .. code-block:: console

          $ aiidalab install quantum-espresso

       This will install the most recent version of an app, regardless of whether it is already installed or not. You will be prompted to confirm the operation.

       ---

       You can install a specific version by using standard `PEP 440 version specifiers`_, for example:

       .. code-block:: console

          $ aiidalab install quantum-espresso==v22.01.0

.. _PEP 440 version specifiers: https://www.python.org/dev/peps/pep-0440/#version-specifiers
