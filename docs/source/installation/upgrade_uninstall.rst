.. _upgrade:

Upgrade or uninstall the app
============================

Upgrading or uninstalling the Quantum ESPRESSO app can be achieved via the App Manager tool or the built-in terminal.

.. tab-set::

   .. tab-item:: App Manager

      .. grid:: 1
         :gutter: 3

         .. grid-item-card:: **Step 1: Open the app management page**

            Click on the **Manage App** button to open the app manager.

            .. image:: ../_static/images/app-management-start-page-upgrade-available.png

         .. grid-item-card:: **Step 2: Upgrade or uninstall the app**

            A green :fa:`arrow-circle-up` **Update** button indicates that a newer version is available. Click it to upgrade the app. Alternatively, you may choose to uninstall the app by clicking on the red :fa:`trash` **Uninstall** button.

            .. image:: ../_static/images/app-management-upgrade-available.png

            .. note::

               By default, the app will be upgraded to the latest available version. Alternatively, you can select any available version including versions lower than the currently installed one.

   .. tab-item:: Terminal

      .. grid:: 1
         :gutter: 3

         .. grid-item-card:: You can upgrade the app from the built-in terminal with

            .. code-block:: console

               $ aiidalab install quantum-espresso

            This will install the most recent version of an app, regardless of whether it is already installed or not. You will be prompted to confirm the operation.

         .. grid-item-card:: You can install a specific version by using standard `PEP 440 version specifiers`_, for example:

            .. code-block:: console

               $ aiidalab install quantum-espresso==v22.01.0

         .. grid-item-card:: You can uninstall the app from the built-in terminal with

            .. code-block:: console

               $ aiidalab uninstall quantum-espresso

            Again, you will be prompted to confirm the operation.

.. _PEP 440 version specifiers: https://www.python.org/dev/peps/pep-0440/#version-specifiers
