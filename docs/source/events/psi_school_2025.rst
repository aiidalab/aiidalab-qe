PSI AiiDAlab Course (7–9 April 2025)
====================================

Overview
--------

This short course at PSI (7–9 April) will introduce participants to
density-functional theory (DFT) simulations using the AiiDAlab Quantum
ESPRESSO app. Topics include fundamental concepts of electronic structure
calculations such as k-point sampling, energy cutoffs, and pseudopotentials.

During the hands-on sessions, participants will carry out structure
relaxations, calculate projected density of states (PDOS), analyze
electronic band structures, and explore advanced materials properties using
various AiiDAlab plugins (e.g., muon spectroscopy, phonons with IR/Raman/INS
spectroscopy, and Bader charge analysis).

All sessions will be live-streamed for remote participants (listen-only),
who are encouraged to follow along with the demonstrations. Links to relevant
AiiDAlab resources will be provided.

Lecturers
---------

- **Prof. Dr. Nicola Marzari** (EPFL, PSI)
- **Dr. Giovanni Pizzi** (PSI)

Schedule and Content
--------------------

**Day 1: Monday, 7 April, 14:00–16:00 (PSI Auditorium)**

- **14:00–15:00**: Basics of Electronic Structure Calculations

  - Overview of k-points and cutoffs
  - Introduction to the SSSP pseudopotential library
  - Demonstration of `seekpath`, `osscar`, and other web tools

- **15:00–16:00**: AiiDAlab Demo & Exercises

  - Logging in to the AiiDAlab demo server, quick demonstration of QE workflows
  - Using ``aiidalab-launch`` for local installations
  - Running crystal relaxations, PDOS, and fat-band calculations
  - Hands-on exercises (guided by in-app tutorials and short videos)
  - Discussion of computational cost (resources, HPC, etc.)

**Day 2: Wednesday, 9 April, 14:00–16:00 (PSI Auditorium)**

- **14:00–14:15**: AiiDAlab Ecosystem Overview

  - Capabilities and plugins
  - Setting up additional HPC resources (Merlin, local clusters, etc.)

- **14:15–15:45**: Advanced Exercises (choose your own)

  - **Phonons + IR/Raman**
  - **Phonons + INS**
  - **aiida-muon**
  - **Bader charges**

- **15:45–16:00**: Outlook and Next Steps

  - Additional features (e.g., XPS workflows)
  - Long-term support and how to continue using AiiDAlab at PSI

Getting Started
---------------

1. **Prepare Your Computer**

   - For hands-on exercises, please bring a laptop with a modern browser (Chrome, Firefox, or Safari).
   - Ensure you have a stable internet connection (on-site or remote).

2. **Access AiiDAlab**

   - **Option 1 (Cloud Access)**: Use the AiiDAlab demo server. Please register and log in here: `AiiDAlab demo`_
   - **Option 2 (Local Installation)**: Use the ``aiidalab-launch`` tool by following the instructions provided here: `AiiDAlab-QE local installation`_

3. **Reading Material**

   - Familiarize yourself with the `AiiDAlab-QE tutorials`_ (we recommendstarting with the basic relaxation workflow).
   - Explore advanced materials property calculations in the `AiiDAlab-QE how-to`_ sections.
   - Review the `SSSP website`_ for details on pseudopotentials.
   - Check out `seekpath`_ or `osscar`_ for geometry and band-structure exploration tools.

4. **Remote Participation**

   - Sessions will be streamed for remote participants (listen-only).
   - You will be able to follow the presentations and watch demos live.

.. _AiiDAlab-QE tutorials: https://aiidalab-qe.readthedocs.io/tutorials/index.html
.. _AiiDAlab-QE how-to: https://aiidalab-qe.readthedocs.io/howto/index.html
.. _AiiDAlab-QE local installation: http://127.0.0.1:8000/installation/access_aiidalab/container.html#aiidalab-launch
.. _AiiDAlab demo: https://demo.aiidalab.io/
.. _SSSP website: https://www.materialscloud.org/discover/sssp
.. _seekpath: https://seekpath.readthedocs.io
.. _osscar: https://osscar.readthedocs.io

Contact and Support
-------------------

If you encounter any issues with your setup or have questions regarding
the tutorials, please contact:

- **Dr. Miki Bonacci** (miki.bonacci@psi.ch)
- **Dr. Xing Wang** (xing.wang@psi.ch)
- **Mr. Timo Reents** (timo.reents@psi.ch)

We look forward to seeing you at the course!
