==========================================================================
 BattMo is a framework for continuum modelling of electrochemical devices.
==========================================================================

.. image:: https://zenodo.org/badge/410005581.svg
   :target: https://zenodo.org/badge/latestdoi/410005581
.. image:: https://github.com/BattMoTeam/BattMo/actions/workflows/runGitHubTests.yml/badge.svg
.. image:: https://github.com/BattMoTeam/BattMo/actions/workflows/doc.yml/badge.svg

The Battery Modelling Toolbox (**BattMo**) is a resource for continuum modelling of electrochemical devices in MATLAB. The initial development features a pseudo X-dimensional (PXD) framework for the Doyle-Fuller-Newman model of lithium-ion battery cells. However, the development plan for BattMo includes extensions to other battery chemistries (e.g. metal-air) and eventually hydrogen systems (i.e. electrolyzers and fuel cells).

**BattMo** offers users a flexible framework for building fully coupled electrochemical-thermal simulations of electrochemical devices using 1D, 2D, or 3D geometries. **BattMo** is implemented in MATLAB and builds on the open-source MATLAB Reservoir Simulation Toolbox (`MRST <https://www.sintef.no/Projectweb/MRST/>`_) developed at SINTEF. MRST provides a solid basis for finite volume grid generation of complex geometries and advanced numerical solvers that enable fast simulations for large systems.

For the latest information including video tutorials and project gallery, please visit the project webpage:
`https://batterymodel.com <https://batterymodel.com/>`_

The documentation is found at the `documentation webpage <https://battmoteam.github.io/BattMo/>`_. We try to do our best to keep it up-to-date.

.. raw:: html

   <img src="Documentation/battmologo_text.png" style="margin-left: 5cm" width="300px">

Installation
------------

Before cloning this reposity you must make sure you have **Git LFS** installed. See `https://git-lfs.com` for instructions on downloading and installation.

BattMo is based on `MRST <https://www.sintef.no/Projectweb/MRST/>`_, which provides a general unstructured grid format,
generic MATLAB automatic differentiation tools and Newton solvers. The MRST source code wil be installed directly via
**git submodules**. To install BattMo, you have therefore to clone this repository with the submodule option
``--recurse-submodules``, as follows:

``git clone --recurse-submodules https://github.com/BattMoTeam/BattMo.git``

Then start MATLAB and in the directory where you cloned the repository, run:

``startupBattMo``

You can check that that your installation is setup correctly by running one of the example scripts:

``runBatteryP2D``

Iterative solvers
-----------------

Iterative solvers are needed to solve large problems with many degrees
of freedom. The 2012 **open source** version of the `AGMG
<http://agmg.eu/>`_ iterative solver is provided as a
`submodule <https://github.com/BattMoTeam/agmg>`_, as well as `AMGCL
<https://github.com/ddemidov/amgcl>`_ are supported.

Tutorials
---------

Tutorials are presented in `documentation <https://BattMoTeam.github.io/BattMo/>`_.

Acknowledgements
-----------------
BattMo has received funding from the European Union’s Horizon 2020 innovation program under grant agreement numbers:

* 875527 HYDRA
* 957189 BIG-MAP
* 101104013 BATMAX
* 101103997 DigiBatt
