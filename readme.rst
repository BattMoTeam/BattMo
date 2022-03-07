========================================================================
 BatMo is framework for continuum modelling of electrochemical devices.
========================================================================

The Battery Modelling Toolbox (**BatMo**) is a resource for continuum modelling of electrochemical devices in MATLAB. The initial development features a pseudo X-dimensional (PXD) framework for the Doyle-Fuller-Newman model of lithium-ion battery cells. However, the development plan for BatMo includes extensions to other battery chemistries (e.g. metal-air) and eventually hydrogen systems (i.e. electrolyzers and fuel cells).

**BatMo** offers users a flexible framework for building fully coupled electrochemical-thermal simulations of electrochemical devices using 1D, 2D, or 3D geometries. **BatMo** is implemented in MATLAB and builds on the open-source MATLAB Reservoir Simulation Toolbox (MRST) developed at SINTEF. MRST provides a solid basis for finite volume mesh generation of complex geometries and advanced numerical solvers that enable fast simulations for large systems.

For the latest information including video tutorials and project gallery, please visit the project webpage:  `https://batterymodel.com <https://batterymodel.com/>`_

.. raw:: html

   <img src="Documentation/batmologo_text.png" style="margin-left: 5cm" width="300px">

Installation
------------

BatMo is based on `MRST <https://www.sintef.no/Projectweb/MRST/>`_, which provides a general unstructured grid format,
generic MATLAB automatic differentiation tools and Newton solvers. The MRST code source wil be installed directly via
**git submodules**. To install batmo, you have therefore to clone this repository with the submodule option
``--recurse-submodules``, as follows:

``git clone --recurse-submodules  git@github.com:batmoTeam/batmo.git``

Then start MATLAB and in the directory where you cloned the repository, run:

``startup``

You can check that that your installation is setup correctly by running one of the example scripts :

``runBattery1D``

Iterative solvers
-----------------

Iterative solvers are needed to solve large problems with many degrees of freedom. The **open source** version from 2012
of the `AGMG <http://agmg.eu/>`_ iterative solver from **2012** is provided as a `submodule
<https://github.com/batmoTeam/agmg>`_. We plan to integrate newer open source iterative solvers such as `AMGCL
<https://github.com/ddemidov/amgcl>`_

Tutorials
---------

Tutorials are presented in `documentation <https://batmoteam.github.io/batmo-doc/>`_ (in progress ...)

Naming Conventions (TBC)
------------------------
Class names are nouns in UpperCamelCase.  
Instance names are nouns in lowerCamelCase.  
Function names are verbs or phrases in lowerCamelCase.  
Common variable names are represented by Latin letters (case set according to convention) or spelled-out lowercase Greek letters (e.g. phi).  
Other variable names may be nouns in lowerCamelCase.  

Conntributors, in alphabetical order
-----------------------------------

* Dr. Simon Clark, SINTEF Industry  
* Dr. Mike Gerhardt, SINTEF Industry  
* Dr. Halvor Møll Nilsen, SINTEF Digital
* Dr. Xavier Raynaud, SINTEF Digital  
* Dr. Roberto Scipioni, SINTEF Industry  

Acknowledgements
-----------------
BatMo has received funding from the European Union’s Horizon 2020 innovation program under grant agreement numbers:

* 875527 HYDRA  
* 957189 BIG-MAP  
