===============
Access AiiDAlab
===============

.. div:: sd-fs-5

   Please select one of the following options for accessing AiiDAlab.

.. grid:: 1 1 1 2
   :gutter: 3

   .. grid-item-card:: Local Docker Instance
      :text-align: center
      :shadow: md

      Install Docker locally and run an instance of an AiiDAlab image *pre-configured** for the Quantum ESPRESSO app. No prior knowledge of Docker necessary!

      ++++

      .. button-ref:: docker
         :ref-type: doc
         :click-parent:
         :expand:
         :color: primary
         :outline:

         To the guide

   .. grid-item-card:: Virtual Machine Image
      :text-align: center
      :shadow: md

      Download a virtual machine image for AiiDAlab based on Quantum Mobile, *pre-configured** with everything you need to run the Quantum ESPRESSO app.

      ++++

      .. button-ref:: vm
         :ref-type: doc
         :click-parent:
         :expand:
         :color: primary
         :outline:

         To the download page

   .. grid-item-card:: Materials Cloud AiiDAlab Server
      :text-align: center
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
      :text-align: center
      :shadow: md

      For members of the Materials Modeling MarketPlace, log into the open AiiDAlab server hosted on the Materials MarketPlace.

      ++++

      .. button-link:: https://aiidalab.materials-marketplace.eu/hub/login
         :click-parent:
         :expand:
         :color: primary
         :outline:

         Launch the server

.. div::

\* The Quantum ESPRESSO app is included in these images. No further installation required!

.. toctree::
   :maxdepth: 1
   :hidden:

   docker
   launch
   vm
