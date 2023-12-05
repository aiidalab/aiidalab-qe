===============
AiiDAlab launch
===============

.. note::

   If you ran the container directly using the :code:`docker` command, feel free to skip this page.

.. important::

   The following steps require a local installation of Docker. If you have yet to do so, you can follow the instructions :doc:`here <docker>`.

`AiiDAlab launch`_ is a thin Docker wrapper which takes care of all the prerequisites to run the AiiDAlab Docker image. To use AiiDAlab launch, you will have to

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
