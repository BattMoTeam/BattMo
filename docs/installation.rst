========================
Installation instruction
========================

BATMO is based on `MRST`_ which provides a general unstructured grid format, generic MATLAB automatic differentiation
tools and Newton solvers.

The MRST code source can be installed directlu using git submodules:

.. code-block:: shell

   git clone --recurse-submodules -b newmaster git@bitbucket.org:mrst/project-batman.git


Then start Matlab and in the directory :code:`project-batman` where you cloned the repo, run:

.. code-block:: matlab

   startup



You can check that that your installation is setup correctly by running one of the example scripts

.. code-block:: matlab

   test_euler

.. _MRST: https://www.sintef.no/Projectweb/MRST/
