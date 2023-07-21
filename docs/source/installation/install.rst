===============
Install the app
===============

.. note::

   This section is mainly for users who accessed AiiDAlab via the open servers (Materials Cloud or Materials MarketPlace), where the app is not yet installed. For all other options, we've taken care of installing the app. You may proceed to use the app. To learn more about how to do so, you can check out the :doc:`tutorials </tutorials/index>`.

   Nevertheless, if interested, please feel free to read on. This guide may be followed in general for other AiiDAlab apps, with minor adjustments if any.

AiiDAlab offers two options for installing the Quantum ESPRESSO app:


#. Via its :fa:`puzzle-piece` **App Store**

#. Via its built-in :fa:`terminal` **Terminal**

.. tip::

   If you're using `aiidalab-launch` to start AiiDAlab locally, you may also install the app directly from the host computer with

   .. code-block::

      aiidalab-launch exec -- aiidalab install quantum-espresso

.. tab-set::

   .. tab-item:: App Store

      .. grid:: 1
         :gutter: 3

         .. grid-item-card:: **Step 1: Open the App Store**

            Click on the :fa:`puzzle-piece` icon in the nav bar. This will open the app store page in a new window or tab.

            .. image:: ../_static/images/nav-bar-app-store.png


         .. grid-item-card:: **Step 2: Search for app**

            Scroll down until you find the app.

            .. tip::

               You can filter apps by category.

            .. image:: ../_static/images/app-management-app-store.png

            Apps not yet installed will show up as follows:

            .. image:: ../_static/images/app-management-app-not-installed.png

         .. grid-item-card:: **Step 3: Install the app**

            Click on the **Install** button to install the app and its dependencies.

            .. note::

               In some cases, you may install a prerelease by clicking on the *Include prereleases* check box. Use this option only if you require access to a not yet released feature, or if you would like to test a new app version and provide feedback to the developer(s).

         .. grid-item-card:: **Step 4: Wait for the installation process to complete**

            The current installation process will be displayed in a terminal widget.

            .. image:: ../_static/images/app-management-app-installation-completed.png

         .. grid-item-card:: **Step 5: Start the app from the start page**

            When the installation process is finished, the newly installed app will show up on the start page. Launch the app by clicking on the Quantum ESPRESSO logo.

            .. image:: ../_static/images/app-management-start-page.png


   .. tab-item:: Terminal

      .. grid:: 1
         :gutter: 3

         .. grid-item-card:: **Step 1: Open the Terminal**

            Open the built-in terminal by clicking on the :fa:`terminal` icon in the nav bar.

            .. image:: ../_static/images/nav-bar-terminal.png

         .. grid-item-card:: **Step 2: Install the app with the aiidalab command line tool**

            .. code-block:: console

               $ aiidalab install quantum-espresso

         .. grid-item-card:: **Step 3: Start the app from the start page**

            When the installation process is finished, the newly installed app will show up on the start page. Launch the app by clicking on the Quantum ESPRESSO logo.

            .. image:: ../_static/images/app-management-start-page.png

.. _AiiDAlab app store: https://aiidalab.github.io/aiidalab-registry
