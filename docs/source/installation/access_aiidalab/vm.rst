=====================
Virtual Machine Image
=====================

.. note::

   If you uses another method to access AiiDAlab, you may proceed to :doc:`install the app </installation/install>`.

The AiiDAlab Quantum ESPRESSO virtual image based on `Quantum Mobile <https://quantum-mobile.readthedocs.io/>`_ is available for the following architectures:

Intel x86_64
------------

Get AiiDAlab Quantum ESPRESSO virtual machine running on your computer in three simple steps:

#. Step1: Download virtual machine image (4.09 GB)

   + URL: https://bit.ly/3Nyf4Dt
   + Filename: `quantum_mobile_23.10.0-qeapp.ova`
   + MD5 hash: `5c80d9ab2458f2edac30e97dc2fe36e7`

#. Step2: Install Virtual Box 6.1.6 or later (see https://www.virtualbox.org)
#. Step3: Import virtual machine image into Virtualbox (11.2 GB) File => Import Appliance

Login credentials:

+ username: `max`
+ password: `moritz`

The default configuration of 4 cores and 4096 MB RAM can be adjusted in the VM settings.

M1/M2 apple silicon
-------------------

If you have an M1/M2 apple silicon, you can use `UTM <https://mac.getutm.app/>`_ to run the AiiDAlab Quantum ESPRESSO virtual machine. The steps are as follows:

#. Step1: Download compressed virtual machine image (5.9 GB)

   + URL: https://bit.ly/486rKtP
   + Filename: `quantum_mobile_23.10.0-qeapp.utm.zip`
   + MD5 hash: `44ea9189d788737459c31bf330366926`

#. Step2: Decompress the zip file.
#. Step3: Install UTM (see https://mac.getutm.app/)
#. Step4 Import image into UTM by clicking on the Create a New Virtual Machine button and selecting the `quantum_mobile_23.10.0-qeapp.utm` file.

Login credentials:

+ username: `max`
+ password: `moritz`

The default configuration of 4 cores and 4096 MB RAM can be adjusted in the VM settings.
