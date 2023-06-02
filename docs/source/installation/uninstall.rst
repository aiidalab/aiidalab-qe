.. _uninstall:

=================
Uninstall the app
=================

Uninstalling an app works similar to upgrading or downgrading an app via the **Manage App** page or on the Terminal.

.. tabbed:: App Manager

    .. panels::
       :container: container-lg pb-3
       :column: col-lg-12 p-2

       **Step 1: Find the app you would like to uninstall on the start page.**

       On the home app start page, simply look for the app you would like to uninstall.

       .. image:: ../_static/images/app-management-start-page.png

       Click on the **Manage App** button to open the app manager.

       ---

       **Step 2: Uninstall**

       The app manager allows you to uninstall the app or to install a different version.

       .. image:: ../_static/images/app-management-app-installed.png

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
