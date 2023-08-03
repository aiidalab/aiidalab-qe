==============
Install Docker
==============

.. note::

   If you used another method to access AiiDAlab, you may proceed to :doc:`install the app </installation/install>`.

AiiDAlab is available as a Docker container - a self-contained, pre-configured environment including all the necessary software to access the AiiDAlab platform. To run the container, you first need to install Docker on your local machine. If you have yet to do so, you may click `here <https://docs.docker.com/get-docker>`_ to follow the official Docker installation guide.

.. important::

   On Linux, you need `root` privileges to perform the `post-installation steps for Docker Engine <https://docs.docker.com/engine/install/linux-postinstall/>`_.

Once Docker is installed, you can launch the container in one of several ways depending on your operating system. This is discussed in the following sections.

Linux
*****

You have two options to launch the container:

#. Execute the following from terminal:

   .. code-block:: console

      docker run -p 8888:8888 aiidalab/qe:edge

#.  Use the `aiidalab-launch <launch>` (**recommended**)

   .. code-block:: console

      pipx install aiidalab-launch # install the aiidalab-launch tool
      aiidalab-launch profile add --image aiidalab/qe:edge QE # add the profile named QE and set to using image aiidalab/qe:edge. Select `n` for the question `Do you want to edit it now? [Y/n]:`
      aiidalab-launch start -p QE

   The Quantum ESPRESSO app will be available at the URL printed in the terminal. Go to the URL in your browser to access the app!

   For more information on the aiidalab-launch tool, please refer to the :doc:`aiidalab-launch documentation <launch>`.

Windows/Mac
***********

.. note::

   Instructions for using the `Docker Desktop <https://www.docker.com/products/docker-desktop>`_ app coming soon
