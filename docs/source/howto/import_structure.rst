.. _import_structure:

================
Import structure
================

You can select a structure from four sources.

- Upload file
- OPTIMADE
- AiiDA database
- From Examples

.. figure:: /_static/images/howto_select_structure.png
   :width: 20cm


Upload
=======

User can upload a structure file directly from the machine. The supported formats are listed in the https://wiki.fysik.dtu.dk/ase/ase/io/io.html#ase.io.read.

From OPTIMADE
=============
`OPTIMADE`_ is a common REST API for accessing structural information of databases that contain calculated properties of existing and hypothetical materials. The OPTIMADE structure importer allows you to search for structures in the a wide range of databases (https://providers.optimade.org/) that support the OPTIMADE API.

In order to select a compound from the OPTIMADE structure importer, you need to

- select the `provider` from the dropdown menu
- select a database for the selected provider (if there are multiple databases available)
- enter the formula of the compound in the search field. Or, select the element from the periodic table.
- click on the `Search` button
- The importer will then search for the compound in the selected database and return a list of structures that match the search criteria. You can then select the structure you want to import from the list.

AiiDA database
==============

Here you can select the structure in your AiiDA database.


From Examples
==============
Here you can select the structure from the examples provided by the `Qeapp` for testing purposes.



.. _OPTIMADE: https://www.optimade.org/
