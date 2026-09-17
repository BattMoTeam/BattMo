============
Installation
============

Release installation
====================

.. _installation:

The latest release version of BattMo is available as a **zip file**. We expect it to work on MATLAB version R2022a or newer.

The latest release version of BattMo is available `here <https://github.com/BattMoTeam/BattMo/releases/latest>`__ as a **zip file**.

1. Create a directory where you want BattMo to be installed. Then, **download** :code:`battmo.zip` in this directory and **unzip** the file

2. Start MATLAB and run the file :code:`startupBattMo` which is located at the root of your BattMo installation directory

   .. code-block:: matlab

      startupBattMo

BattMo is now **installed**. You can check that your installation is setup correctly by running one of the example scripts, directly from Matlab command line.

.. code-block:: matlab

   runBatteryP2D


Installation from git
=====================

To install the development version of BattMo from git, follow the
`installation instructions in the repository README
<https://github.com/BattMoTeam/BattMo/blob/main/readme.rst#installation>`_.
These cover Git LFS setup, cloning with submodules, and starting BattMo in MATLAB.

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


 Contributing
 ============

 To contribute, please see `contributing.rst <https://github.com/BattMoTeam/BattMo/blob/main/contributing.rst#>`_.
