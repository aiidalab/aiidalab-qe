=====================
Local Docker Instance
=====================

.. note::

   If you used another method to access AiiDAlab, you may proceed to :doc:`install the app </installation/install>`.

AiiDAlab is available as a Docker container - a self-contained, pre-configured environment including all the necessary software to access the AiiDAlab platform.
Conveniently, we provide a bluprint (image) for such a container with the Quantum ESPRESSO app pre-installed and ready for use. To run the container, you first need to `install Docker <https://docs.docker.com/get-docker>`_ on your local machine.

.. important::

   The Docker installation link above walks you through installing Docker Desktop - a convenient graphical user interface to Docker. However, if you have chosen instead to install the `docker engine <https://docs.docker.com/engine/install/ubuntu/>`_ directly for use via a terminal, if **(and only if)** you are on a Linux system, you will need `root` privileges to perform the `post-installation steps for Docker Engine <https://docs.docker.com/engine/install/linux-postinstall/>`_.

Once Docker is installed, you can launch the container in one of several ways depending on your operating system. This is discussed in the following sections.

Docker Desktop
**************

Note that Docker Desktop is also available for MacOS and Linux. However, you may prefer to use the command line interface (CLI). If so, proceed to the following sections for instructions.

AiiDAlab launch
***************

.. important::

   The following steps require a local installation of Docker. You can verify your Docker installation by running ``docker run hello-world`` in the terminal.

`AiiDAlab launch`_ is a thin Docker wrapper which takes care of all the prerequisites to run the AiiDAlab Docker image. It helps to manage multiple AiiDAlab profiles, each with its own home directory for persistent storage, and allows to easily switch between them.
To use AiiDAlab launch, you will have to

#. Install AiiDAlab launch with `pipx <https://pypa.github.io/pipx/installation/>`_ (**recommended**):

   .. code-block:: console

      pipx install aiidalab-launch

   or directly with pip

   .. code-block:: console

      pip install aiidalab-launch

#. Set up a new `QE` profile with

   .. code-block:: console

      aiidalab-launch profile add --image aiidalab/qe:latest QE

   At the prompt, enter `n` to skip editing the profile settings.

#. Start AiiDAlab with

   .. code-block:: console

       aiidalab-launch start -p QE

#. Follow the URL on the screen to open AiiDAlab in the browser

.. tip::

   For more detailed help, run

   .. code-block:: console

      aiidalab-launch --help

Profile Management
^^^^^^^^^^^^^^^^^^

As shown above, you can manage multiple profiles in AiiDAlab launch, e.g., with different home directories or ports. For more information, run

.. code-block:: console

   aiidalab-launch profiles --help

You can inspect the status of all configured AiiDAlab profiles with

.. code-block:: console

   aiidalab-launch status

.. _`AiiDAlab launch`: https://github.com/aiidalab/aiidalab-launch

Using docker CLI directly
*************************

It is not necessary to use AiiDAlab launch to run the AiiDAlab container.
You can also use the docker CLI directly by running

.. code-block:: console

   docker run -p 8888:8888 aiidalab/qe:latest

Follow the URL on the screen to open AiiDAlab in the browser.

.. important::

   If you use the docker CLI directly, the data in the home directory of the container will be lost when the container is deleted. You can use the ``-v`` option to mount a local directory to the container to store the data persistently. For more information, run ``docker run --help``.
