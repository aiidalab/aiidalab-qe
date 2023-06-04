.. _how_to_setup_computer_code:

Setup remote computer/code
==========================

The instruction in collapsible text is for CSCS machine.

Now you need set up a new computer (e.g., the CSCS Daint or Eiger machine if you have acount on it) and the codes on this newly setup computer (i.e., you need to tell AiiDA where the executables are, which modules to load etc.).

.. note::

   The widget to assist you in the computer/code setup is not very convinient for CSCS machines, because CSCS requires to specify the account and setting up a default one is not currently straightforward in AiiDA; we are working to improve this.


**Setup new code**
- Click the button `Setup new code`.
- In the `Domain` field, please select `daint.cscs.ch`.
- Click `Detailed Setup`, then it will show three steps:

    - password-less SSH
    - computer in AiiDA
    - code in AiiDA

**Set up password-less SSH connection**
- In the `SSH username`, input your username for CSCS.
- Click `Setup ssh`.
- On the textbox above `Detailed Setup`, input your password, and click `Continue`.
- Wait a few seconds, it will ask you to input your password again, and click `Continue`

If you see a message as below:
```
30s timeout ...
```
Please click the `Setup ssh` button again and input the password.


**Setup a computer in AiiDA**

Click the button `Setup computer`, then add one line in the `Prepend text` to specify your account (it is the project account, such as `mr32` or `mrcloud`) at CSCS. Replace the `<account_name>` string with it!
```
#SBATCH --account=<account_name>
```
You can find further parameters in the screenshot below.

.. figure:: /_static/images/setup_computer.png
   :width: 12cm


Click the button `test computer`, then wait 1 minute: you will see this log message if everything was setup correctly:

.. code-block:: none

   Report: Testing computer<daint-mc> for user<aiida@localhost>...
   * Opening connection... [OK]
   * Checking for spurious output... [OK]
   * Getting number of jobs from scheduler... [OK]: 2255 jobs found in the queue
   * Determining remote user name... [OK]: xingwang
   * Creating and deleting temporary file... [OK]
   Success: all 5 tests succeeded


**Set up a code**

In the field `Select computer`, please select the computer you just setup, in this case `daint-mc`. Check all other options, and then click the button `Setup code` to setup the `pw.x` code.

.. figure:: /_static/images/setup_code.png
   :width: 12cm

Similarly, setup the `dos.x` and `projwfc.x` codes, needed by the app to compute DOS and PDOS. You can skip the computer setup step, since it needs to be done only once per comouter, and go to `Set up a code in AiiDA` directly to setup the remaining codes.

After you finishing the codes setup, you can launch a new calculation with the newly setup codes, that should be called:

- pw-7.0@daint-mc
- dos-7.0@daint-mc
- projwfc-7.0@daint-mc

.. figure:: /_static/images/select_new_code.png
   :width: 12cm
