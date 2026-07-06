============
Installation
============

Release installation
====================

.. _installation:

The latest release version of BattMo is available `here <https://github.com/BattMoTeam/BattMo-test/releases/latest>`__ as a **zip file**.

1. **Download** :code:`battmo.zip` in the directory of your choice and **unzip** the file there.

2. Start MATLAB and run the file :code:`startupBattMo` which is located at the root of the directory

   .. code-block:: matlab

      startupBattMo

BattMo is now **installed**. You can check that that your installation is setup correctly by running one of the example scripts, directly from Matlab command line.

.. code-block:: matlab

   runBatteryP2D

   
Installation from git
=====================

BattMo source code can be installed using git. In this way, you can easily keep track of the last developments.

First, open a terminal and navigate to the directory where you would like to install BattMo. Then, clone the repository using the following command, which will include all the dependencies as submodules

.. code-block:: shell

   git clone --recurse-submodules https://github.com/BattMoTeam/BattMo.git

Then, run :code:`startupBattMo` 

Here is video which guides you through the installation in details

.. youtube:: -XVppzyNSs0

For a detailed guided installation of git, you can consult this `video <https://www.youtube.com/watch?v=FMXpZjXhaFY>`_

.. _MRST: https://www.sintef.no/Projectweb/MRST/


Update existing installation
============================

In the case where we alread have installed BattMo and you want to update to the latest version. As usual in git, you
will do that by running

.. code-block:: shell

   git pull

In addition to that, the dependencies that are given through git submodules. The are not updated often but, if it is the
case, you will need to run in addition to the previous command,

.. code-block:: shell

   git submodule update --recursive

