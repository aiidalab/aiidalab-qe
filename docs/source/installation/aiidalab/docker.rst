=================
Installing Docker
=================

.. note::

   If you used another method to access AiiDAlab, feel free to skip this page.

AiiDAlab is available as a Docker container - a self-contained, pre-configured environment including all the necessary software to access the AiiDAlab platform. To run the container, you first need to install Docker on your local machine. If you have yet to do so, you may click `here <https://docs.docker.com/get-docker>`_ to follow the official Docker installation guide.

.. important::

   On Linux, you need `root` privileges to perform the `post-installation steps for Docker Engine <https://docs.docker.com/engine/install/linux-postinstall/>`_.

Once Docker is installed, you have two options to launch the container:

#. Execute the following from terminal:

   .. code-block:: console

      docker run -p 8888:8888 aiidalab/aiidalab-docker-stack/qe

#.  Use the :doc:`aiidalab-launch <aiidalab_launch>` tool (**recommended**)
