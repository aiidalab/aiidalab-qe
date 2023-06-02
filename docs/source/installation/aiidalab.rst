.. _access_aiidalab:

===============
Access AiiDAlab
===============

There are three options for accessing AiiDAlab:

1. Log into one of the `open AiiDAlab servers <https://www.aiidalab.net/deployments/>`_
2. Run ``aiidalab`` from the `Quantum Mobile Virtual Machine <https://quantum-mobile.readthedocs.io/>`_ terminal
3. Run the AiiDAlab docker container locally using :ref:`aiidalab-launch <aiidalab_launch>` (see below)

.. _aiidalab_launch:

***************
AiiDAlab launch
***************

To run AiiDAlab on your own workstation or laptop you can either

 - (**recommended**) Use the `aiidalab-launch <https://github.com/aiidalab/aiidalab-launch#aiidalab-launch>`_ tool which is a thin docker wrapper, or
 - run the image directly with: ``docker run -p 8888:8888 aiidalab/full-stack``.

To use **AiiDAlab launch** you will have to

#. `Install Docker on your workstation or laptop. <https://docs.docker.com/get-docker/>`_

   .. note::

      If you are using Linux, you need to have `root` privileges to do `post-installation steps for the Docker Engine <https://docs.docker.com/engine/install/linux-postinstall/>`_.

#. Install AiiDAlab launch with `pipx <https://pypa.github.io/pipx/installation/>`_ (**recommended**):

   .. code-block:: console

      pipx install aiidalab-launch

   Or directly with pip (``pip install aiidalab-launch``).

#. Start AiiDAlab with

   .. code-block:: console

       aiidalab-launch start

#. Follow the instructions on screen to open AiiDAlab in the browser.

See ``aiidalab-launch --help`` for detailed help.

Instance Management
^^^^^^^^^^^^^^^^^^^

You can inspect the status of all configured AiiDAlab profiles with:

.. code-block:: console

   aiidalab-launch status

Profile Management
^^^^^^^^^^^^^^^^^^

The tool allows to manage multiple profiles, e.g., with different home directories or ports.
See ``aiidalab-launch profiles --help`` for more information.
