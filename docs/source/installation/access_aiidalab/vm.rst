=====================
Virtual Machine Image
=====================

.. note::

   If you used another method to access AiiDAlab, you may proceed to :doc:`install the app </installation/install>`.

The AiiDAlab Quantum ESPRESSO virtual image based on `Quantum Mobile <https://quantum-mobile.readthedocs.io/>`_ image and is available for the following architectures:

+ Intel x86_64
+ Mac computer with Apple Silicon

See the following sections for instructions on how to download and run the virtual machine image on your computer.

Intel x86_64
------------

Get AiiDAlab Quantum ESPRESSO virtual machine running on your computer in three simple steps:

#. Download the virtual machine image (4.09 GB)

   + URL: https://bit.ly/48qUay9 (Google drive: https://bit.ly/3Nyf4Dt)
   + Filename: `quantum_mobile_23.10.0-qeapp.ova`
   + MD5 hash: `5c80d9ab2458f2edac30e97dc2fe36e7`

#. Install Virtual Box 6.1.6 or later (see https://www.virtualbox.org)
#. Import the virtual machine image into Virtualbox (11.2 GB) File => **Import Appliance**

Login credentials:

+ username: `max`
+ password: `moritz`

The default configuration of 4 cores and 4096 MB RAM can be adjusted in the VM settings.

Mac computer with Apple Silicon
-------------------------------

If you have an Apple silicon Mac computer with M1/M2 chips, you can use `UTM <https://mac.getutm.app/>`_ to run the AiiDAlab Quantum ESPRESSO virtual machine. The steps are as follows:

#. Download the compressed virtual machine image (5.9 GB)

   + URL: https://bit.ly/477mMeR (Google Drive: https://bit.ly/486rKtP)
   + Filename: `quantum_mobile_23.10.0-qeapp.utm.zip`
   + MD5 hash: `44ea9189d788737459c31bf330366926`

#. Decompress the zip file.
#. Install UTM (see https://mac.getutm.app/)
#. Import the image into UTM by clicking on the **Create a New Virtual Machine** button and selecting the `quantum_mobile_23.10.0-qeapp.utm` file.

Login credentials:

+ username: `max`
+ password: `moritz`

The default configuration of 4 cores and 4096 MB RAM can be adjusted in the VM settings.
