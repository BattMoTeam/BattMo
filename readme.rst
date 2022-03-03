BatMo is framework for continuum modelling of electrochemical devices.
================================================================

.. raw:: html

   <img src="Documentation/batmologo.png" style="margin-left: 5cm" width="300px">

Installation
------------
BATMO is based on MRST which provides a general unstructured grid format, generic MATLAB automatic differentiation tools and Newton solvers. The MRST code source can be installed directly using git submodules:

``git clone --recurse-submodules  git@github.com:batmoTeam/batmo.git``

Then start MATLAB and in the directory project-batman where you cloned the repo, run:

``startup``

You can check that that your installation is setup correctly by running one of the example scripts

``runBattery1D``

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
|* 875527 HYDRA  
|* 957189 BIG-MAP  
