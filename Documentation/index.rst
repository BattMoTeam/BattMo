.. BattMo documentation master file, created by
   sphinx-quickstart on Thu Sep 30 20:16:46 2021.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. toctree::
   :hidden:

   functioninterface
   ..
      Installation <installation>
      basicusage
      tutorials
      intermediate
      modeling
      architecture
      json
      geometryinput
      juliabridge
      Electrolyser simulation <publishedExamples/runElectrolyser>
      protonicmembrane
      Computational Graph <computationalGraph/graphdoc>
      app
   ..
      seealso
      References <bibliography>


Welcome
=======

Welcome to the Battery Modelling Toolbox (**BattMo**), a comprehensive solution for continuum modelling of electrochemical devices in `MATLAB <https://se.mathworks.com/products/matlab.html>`_ and `Julia <https://julialang.org/>`_!

**BattMo** facilitates a deep understanding of these devices by simulating cell-level performance in a virtual space. It allows you to calculate dynamic spatial profiles for essential quantities like concentration, electric potential, and temperature. Initially, **BattMo** focuses on the Doyle-Fuller-Newman model for lithium-ion battery cells but has a broader development plan that includes extensions to other battery chemistries such as Na-ion, solid-state, metal-air, and zinc-based systems, along with hydrogen systems like electrolyzers and fuel cells.

Our toolbox offers a flexible framework for building fully coupled electrochemical-thermal simulations with 1D, 2D, or 3D geometries. Powered by the open-source MATLAB Reservoir Simulation Toolbox (MRST), **BattMo** provides efficient finite volume grid generation and advanced numerical solvers, ensuring swift simulations even for complex systems. Whether you're a researcher or developer, **BattMo** is your gateway to unlocking the potential of continuum modelling for electrochemical devices. Dive into our documentation and explore the possibilities!

For the latest information including video tutorials and project gallery, please visit the project webpage :
`https://batterymodel.com <https://batterymodel.com/>`_

.. image:: battmologo.png
   :width: 50%
   :align: center
   :target: https://batterymodel.com/

.. note::
  This project is under active development.


Acknowledgements
================

BattMo has received funding from the European Union’s Horizon 2020 and Horizon Europe innovation programs under grant agreement numbers:

- 875527 - Hybrid power-energy electrodes for next-generation lithium-ion batteries (HYDRA)
- 957189 - Battery interface genome and materials acceleration platform (BIG-MAP)
- 101069765 - Innovative and Sustainable High Voltage Li-ion Cells for Next Generation (EV) Batteries (IntelLiGent)
- 101104031 - Battery management by multi-X (X=scale/physics/use/domain) digital twins (BATMAX)
- 101103997 - Digital Solutions for Accelerated Battery Testing (DigiBatt)
