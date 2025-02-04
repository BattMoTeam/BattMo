---
title: 'BattMo'
tags:
  - Battery
  - Numerical simulation
authors:
  - name: Simon Clark
    orcid: 0000-0002-8758-6109
    affiliation: 1
  - name: Xavier Raynaud
    orcid: 0000-0002-4100-3035
    affiliation: 2
  - name: Halvor Møll Nilsen
    orcid: 0000-0002-2153-0962
    affiliation: 2
  - name: August Johansson
    orcid: 0000-0001-6950-6016
    affiliation: 2
  - name: Eibar Flores
    orcid: 0000-0003-2954-1233
  - name: Lorena Hendrix
    orcid: 0009-0006-9621-6122
    affiliation: 1
  - name: Francesca Watson
    orcid: 0000-0002-4391-4166
    affiliation: 2
  - name: Sridevi Krishnamurthi
    orcid: 0009-0006-0805-6713
    affiliation: 1
  - name: Olav Møyner
    orcid: 0000-0001-9993-3879
    affiliation: 2
affiliations:
 - name: SINTEF Industry, Dept. of Sustainable Energy Technology, Norway
   index: 1
 - name: SINTEF Digital, Dept. of Mathematics and Cybernetics, Norway
   index: 2
date: 4 February 2025
bibliography: paper.bib
---

# Summary

New high-performance battery designs are essential to achieve the goals of the electric energy transition. Across the
industrial landscape, companies are adopting battery systems with ever increasing requirements for performance,
lifetime, and cost. Market demand for batteries is projected to reach 390 GWh/year in 2030, and is growing so quickly
that it may strain global supply chains. At the same time, research institutions are generating an abundance of data and
new modelling approaches in the pursuit of novel battery materials and chemistries. Developing rigorous digital
workflows can help industrial and research institutions reduce the need for costly physical prototyping and derive
greater insight and knowledge from their data.

Flexible and adaptable modelling frameworks are a connerstone of battery digitalization. High-quality battery models can
be labor-intensive to develop from scratch. Recently, a variety of open-source battery modelling codes have been
released including PyBaMM [@sulzer2021python], cideMOD [@CiriaAylagas2022], LIONSIMBA [@torchio2016lionsimba], PETLion
[@Berliner_2021], and MPET, among others. These open-source modelling frameworks help the battery community reduce the
cost of model development and help ensure the validity and the reproducibility of findings.

Most physics-based continuum models implement some version of the Doyle-Fuller-Newman approach for Li-ion batteries. The
cell is typically approximated as a one-dimensional mesh with additional discretization along the radius of the active
particles. This pseudo-two-dimensional (P2D) approach has been widely used in the battery modelling community as a
relatively inexpensive way to get greater insight about the evolution of key quantities like concentration, temperature,
and electric potential during cell operation. However, the simplified P2D mesh cannot resolve the effects of the cell
geometry (e.g. current collector tabs). Furthermore, these approaches focus almost exclusively on lithium-ion cell
chemistry.

There is a clear need for a battery modelling framework which can be adapted for both Li-ion and post-Li-ion
technologies and can simulate the performance of battery cells in full three-dimensional designs. Some initial steps in
this direction have been taken. PyBaMM includes a lead-acid battery extension for a python-based P2D framework, while
cideMOD utilizes the Fenics finite element package to simulate Li-ion pouch cells in 3D.

This paper presents the Battery Modelling Toolbox (BattMo), a flexible finite volume continuum modelling framework for
simulating the performance of electrochemical cells. BattMo can quickly setup and solve models for a variety of battery
chemistries, even considering complex designs like cylindical and prismatic cell jelly rolls. Furthermore, the parameter
and result files are supplemented by annotated metadata using the Battery Interface Ontology (BattINFO) to support
semantic interoperability in accordance with the FAIR principles. BattMo builds on the MATLAB Reservoir Simulation
Toolbox (MRST) which provides a reliable foundation for meshing intricate geometries, efficiently solving large systems
of equations, and visualizing the results.


# Statement of need

- Standard data input (json based)
- Different battery format (geometry)
- Flexible tool 

# Physical models

- PXD model (Newman)
- Sea Water model
- Electrolysis models

# Numerical method

Implementation is based on MRST

- Finite volume method
- Implicit time integrator for robustness (work in progress for matlab ode solver)
- Assembly based on automatic differentiation
- Setup for adjoint computation
- AGMG preconditioner and linear solver
- Computational graph model development

# Geometry

![Battery geometries \label{fig:geometries}](figs/batterygeometries.png){width=100%}

# Optimization

# Citations

# Figures

Figures can be included like this:
<!-- ![Caption for example figure.\label{fig:example}](figure.png) -->
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
<!-- ![Caption for example figure.](figure.png){ width=20% } -->

# Acknowledgements

We acknowledge contributions from the EU, Grant agreements 101069765, 875527, 101104013, 101103997

# References
