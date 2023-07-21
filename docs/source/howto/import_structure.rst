.. _import_structure:

================
Import structure
================

You can select a structure from the following four sources:

- Upload file from your computer
- Find a structure using the `OPTIMADE`_ API
- Load an existing structure from the underlying AiiDA database
- Load one structure from a selection of examples

.. figure:: /_static/images/howto_select_structure.png
   :width: 20cm

Upload
======

User can upload a structure file directly from their machine. The supported formats are listed in the https://wiki.fysik.dtu.dk/ase/ase/io/io.html#ase.io.read.

From OPTIMADE
=============

`OPTIMADE`_ is a common REST API for accessing structural information of databases that contain calculated properties of existing and hypothetical materials. The OPTIMADE structure importer allows you to search for structures from the wide range of databases (https://providers.optimade.org/) that support the OPTIMADE API.

In order to select a compound from the OPTIMADE structure importer, you need to

- Select the `provider` from the dropdown menu
- select a (sub)database for the selected provider (if there are multiple databases available)
- Enter the formula of the compound in the search field, or select the element from the periodic table
- Click on the `Search` button
- The importer will then search for the compound in the selected database and return a list of structures that match the search criteria. You can then select the structure you want to import from the list

AiiDA database
==============

`AiiDA`_, the workflow management system serving as the backend to the AiiDAlab platform, will track the provenance of your submitted calculations, including input (output) structures when provided (produced). You can select these structures here from your AiiDA database. This is for instance useful to, e.g., start a calculation from a previously relaxed one, or if you have already imported some structures in AiiDA.

From Examples
=============

Here you can select the structure from examples included with the Quantum ESPRESSO app for testing purposes.

.. _OPTIMADE: https://www.optimade.org/
