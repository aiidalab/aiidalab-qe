.. note::

   If you used another method to access AiiDAlab, you may skip this page.

======
Docker
======

AiiDAlab is available as a Docker container - a self-contained pre-configured environment including all the necessary software. To run the container, you must first install `Docker`_ on your workstation or laptop.

.. important::

   On Linux, you need to have `root` privileges to do the `post-installation steps for Docker Engine <https://docs.docker.com/engine/install/linux-postinstall/>`_.

Once Docker is installed, you have two options to launch the container:

#. Execute the following from terminal:

   .. code-block:: console

      docker run -p 8888:8888 aiidalab/aiidalab-docker-stack/qe

#.  Use the :ref:`aiidalab-launch <aiidalab_launch>` tool (**recommended**)


.. _aiidalab_launch:

AiiDAlab launch
***************

.. important::

   The following steps require a local installation of `Docker`_.

AiiDAlab launch is a thin Docker wrapper which takes care of all the prerequisites to run the AiiDAlab Docker image. To use AiiDAlab launch you will have to

#. Install AiiDAlab launch with `pipx <https://pypa.github.io/pipx/installation/>`_ (**recommended**):

   .. code-block:: console

      pipx install aiidalab-launch

   or directly with pip

   .. code-block:: console

      pip install aiidalab-launch

#. Start AiiDAlab with

   .. code-block:: console

       aiidalab-launch start

#. Follow the instructions on the screen to open AiiDAlab in the browser

.. tip::

   For more detailed help, run

   .. code-block:: console

      aiidalab-launch --help

Profile Management
^^^^^^^^^^^^^^^^^^

You can manage multiple profiles in AiiDAlab launch, e.g., with different home directories or ports. For more information, run

.. code-block:: console

   aiidalab-launch profiles --help

You can inspect the status of all configured AiiDAlab profiles with

.. code-block:: console

   aiidalab-launch status
