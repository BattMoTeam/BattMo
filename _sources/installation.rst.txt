=========================================
Installation and First Steps
=========================================

.. _installation:

Install BattMo by cloning the git repository using :code:`--recurse-submodules` to also install the `MRST`_ dependencies.

First, open a terminal and navigate to the directory where you would like to install BattMo. Then, clone the repository using the following command:

.. code-block:: shell

   git clone --recurse-submodules https://github.com/BattMoTeam/BattMo.git


Start MATLAB and in the directory :code:`battmo` where you cloned the repository. In the MATLAB Command Window run:

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
