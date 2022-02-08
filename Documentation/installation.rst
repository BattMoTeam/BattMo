========================
Installation instruction
========================

BATMO is based on `MRST`_ which provides a general unstructured grid format, generic MATLAB automatic differentiation
tools and Newton solvers.

The MRST code source can be installed directlu using git submodules:

.. code-block:: shell

   git clone --recurse-submodules  git@github.com:batmoTeam/batmo.git


Then start Matlab and in the directory :code:`project-batman` where you cloned the repo, run:

.. code-block:: matlab

   startup


You can check that that your installation is setup correctly by running one of the example scripts

.. code-block:: matlab

   runBattery1D

.. _MRST: https://www.sintef.no/Projectweb/MRST/
