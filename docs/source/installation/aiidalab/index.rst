===============
Access AiiDAlab
===============

.. div:: sd-fs-5

   Please select one of the following options for accessing AiiDAlab.

.. grid:: 1 1 1 2
   :gutter: 3

   .. grid-item-card:: Local Docker instance
      :text-align: left
      :shadow: md

      Walk through installing Docker on your local machine, spinning up a container instance of a pre-configured AiiDAlab image, and running AiiDAlab in your browser. No prior knowledge of Docker necessary!

      ++++

      .. button-ref:: docker
         :ref-type: doc
         :click-parent:
         :expand:
         :color: primary
         :outline:

         To the guide

   .. grid-item-card:: Quantum Mobile Virtual Machine
      :text-align: left
      :shadow: md

      We've pre-configured a Quantum Mobile virtual machine image for AiiDAlab, with everything you need to run the Quantum ESPRESSO app. Download the image here and run it with your favorite virtual machine software.

      ++++

      .. button-link:: https://quantum-mobile.readthedocs.io/
         :click-parent:
         :expand:
         :color: primary
         :outline:

         Download the virtual image

   .. grid-item-card:: Materials Cloud AiiDAlab Server
      :text-align: left
      :shadow: md

      For researchers affiliated with Materials Cloud partners, log into the open AiiDAlab server hosted on the Materials Cloud.

      ++++

      .. button-link:: https://aiidalab.materialscloud.org/hub/login
         :click-parent:
         :expand:
         :color: primary
         :outline:

         Launch the server

   .. grid-item-card:: Materials MarketPlace AiiDAlab Server
      :text-align: left
      :shadow: md

      For members of the Materials Modeling MarketPlace, log into the open AiiDAlab server hosted on the Materials MarketPlace.

      ++++

      .. button-link:: https://materials-marketplace.aiidalab.net/hub/login
         :click-parent:
         :expand:
         :color: primary
         :outline:

         Launch the server

.. toctree::
   :maxdepth: 1
   :hidden:

   docker
   aiidalab_launch
