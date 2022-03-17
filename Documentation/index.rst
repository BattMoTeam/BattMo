.. BattMo documentation master file, created by
   sphinx-quickstart on Thu Sep 30 20:16:46 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Introduction
=================================
The Battery Modelling Toolbox (**Battmo**) is a resource for continuum modelling of electrochemical devices in MATLAB. The initial development features a pseudo X-dimensional (PXD) framework for the Doyle-Fuller-Newman model of lithium-ion battery cells. However, the development plan for BattMo includes extensions to other battery chemistries (e.g. metal-air) and eventually hydrogen systems (i.e. electrolyzers and fuel cells).

**BattMo** offers users a flexible framework for building fully coupled electrochemical-thermal simulations of electrochemical devices using 1D, 2D, or 3D geometries. **BattMo** is implemented in MATLAB and builds on the open-source MATLAB Reservoir Simulation Toolbox (MRST) developed at SINTEF. MRST provides a solid basis for finite volume mesh generation of complex geometries and advanced numerical solvers that enable fast simulations for large systems.

For the latest information including video tutorials and project gallery, please visit the project webpage:  `https://batterymodel.com <https://batterymodel.com/>`_

.. image:: battmologo.png
   :width: 50%
   :align: center
   :target: https://batterymodel.com/

.. note::
  This project is under active development.

BattMo documentation
=================================

.. toctree::
   :maxdepth: 2
              
   Installation <installation>
   tutorials
   inputs
   ..
      batterymodel
   ..
      genmodels
   ..
      electrodes
   ..
      electrolyte
   mrst
   bibliography

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Acknowledgements
=================================
BattMo has received funding from from the European Union’s Horizon 2020 innovation program under grant agreement numbers:

875527 - Hybrid power-energy electrodes for next-generation lithium-ion batteries (HYDRA)
957189 - Battery interface genome and materials acceleration platform (BIG-MAP)
