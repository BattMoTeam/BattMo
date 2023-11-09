=========================================
Installation instruction and requirements
=========================================

.. _installation:

BattMo is based on `MRST`_ which provides a general unstructured grid format, generic MATLAB automatic differentiation
tools and Newton solvers.

The MRST code source can be installed directlu using git submodules:

.. code-block:: shell

   git clone --recurse-submodules https://github.com/BattMoTeam/BattMo.git


Then start Matlab and in the directory :code:`battmo` where you cloned the repo, run:

.. code-block:: matlab

   startupBattMo


You can check that that your installation is setup correctly by running one of the example scripts

.. code-block:: matlab

   runBatteryP2D

Here is video which guides you through the installation in details

.. youtube:: -XVppzyNSs0
              
For a detailed guided installation of git, you can consult this `video <https://www.youtube.com/watch?v=FMXpZjXhaFY>`_
            
.. note::
   
   For large models, iterative solvers are necessary. The **open source** version from 2012 of the `AGMG
   <http://agmg.eu/>`_ iterative solver from **2012** is provided as a `submodule <https://github.com/batmoTeam/agmg>`_. We
   plan to integrate newer open source iterative solvers such as `AMGCL <https://github.com/ddemidov/amgcl>`_

.. _MRST: https://www.sintef.no/Projectweb/MRST/

