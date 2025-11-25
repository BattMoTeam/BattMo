---
title: 'BattMo - Battery Modelling Toolbox'
tags:
  - Battery modeling
  - Numerical simulation
  - MATLAB
authors:
  - name: Xavier Raynaud^[corresponding author]
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
  - name: Simon Clark
	orcid: 0000-0002-8758-6109
	affiliation: 1
affiliations:
 - name: SINTEF Industry, Dept. of Sustainable Energy Technology, Norway
   index: 1
 - name: SINTEF Digital, Dept. of Mathematics and Cybernetics, Norway
   index: 2
date: 24 November 2025
bibliography: paper.bib
---
<!-- To compile this file, after installing docker, from this directory, run : docker run --rm --volume $PWD:/data --user $(id -u):$(id -g) --env JOURNAL=joss openjournals/inara  -->
# Summary

This paper presents the Battery Modelling Toolbox (BattMo), a flexible finite volume continuum modelling framework in MATLAB for
simulating the performance of electrochemical cells. BattMo can quickly setup and solve models for a variety of battery
chemistries, even considering complex designs like cylindical and prismatic cell jelly rolls.

The simulation input parameters, including the material parameters and geometric descriptions, are specified through
JSON schemas. In this respect, we follow the guidelines of the Battery Interface Ontology (BattINFO) to support semantic
interoperability in accordance with the FAIR principles.

The Doyle-Fuller-Newman (DFN) [@Doyle1993ModelingCell] approach is used as a base model. We include fully coupled
thermal simulations. It is possible to include degradation mechanisms such as SEI layer growth) and the use of composite
material, such as a mixture of Silicon and graphite.

The models are setup in a hierarchical way, for clarity and modularity. Each model corresponds to a computational graph,
which indroduces a set of variables (the nodes) and functional relationship (the edges). This design enables the
flexibility for changing and designing new models. This flexibility is also illustrated in the implementation of an
electrolyser model.

The solver in BattMo uses automatic differentiation and support adjoint computation. We can therefore computute the
derivative of objective functions with respect to all parameters in a very efficient way. Gradient-based optimization
routines can be used to calibrate parameters from experimental data by minimizing the difference between observed and
predicted results.

# Statement of need

New high-performance battery designs are essential to achieve the goals of the electric energy transition. Developing
rigorous digital workflows can help industrial and research institutions reduce the need for costly physical prototyping
and derive greater insight and knowledge from their data.

There is a clear need for a battery modelling framework which can be adapted for both Li-ion and post-Li-ion
technologies and can simulate the performance of battery cells in full three-dimensional designs. Recently, a variety of
open-source battery modelling codes have been released including PyBaMM [@sulzer2021python], cideMOD
[@CiriaAylagas2022], LIONSIMBA [@torchio2016lionsimba], PETLion [@Berliner_2021], and MPET, among others. These
open-source modelling frameworks help the battery community reduce the cost of model development and help ensure the
validity and the reproducibility of findings.

BattMo prolongs this effort by supporting fully 3D geometry and includes the possibility to easily modify the underlying
equations. We provide a library of standard battery geometry which is parameterize. Design optimization can also be done
on the geometry, which is an essential part of the design.

A challenge for physically based model for battery is the difficulty to calibrate the parameters. With an adjoint-based
approach, we can effectively calibrate the models from experiments in a reasonable computational time.

# Functionality overview

The default and easiest way to send the input parameters for the simulator is by using a JSON file. The JSON format is a
text based format, which means that the file can be read directly by the user and easily modify. The keywords used in a
BattMo JSON input file are all specified through a set of JSON schema.

List of features

- Standard data input (JSON based with schema)
- Library of battery formats (to be extended)
- Flexible model design
- Visualization
- Parameter calibration
- Design optimization
- Support for standard protocols such as CC, CV, CCCV, and time series
- SEI layer growth model
- Composite material model
- Material database
- Silicon swelling (to be soon merged...)

# Software dependencies

BattMo builds on the MATLAB Reservoir Simulation Toolbox [@mrst-book-i][@MRST:2025] which provides a reliable foundation for meshing
intricate geometries, efficiently solving large systems of equations, and visualizing the results. It is implemented in
MATLAB and seeks to provide Octave compability. Neither BattMo nor MRST rely on extra MATLAB packages; the basic license is sufficient. We recommend using AMG preconditioners from the open-source
AMGCL package [@Demidov2020].

# Battery format library

We support coin cells, jelly roll cells and multipouch cells with different tab layouts. The geometries are parameterized and can be modified
using a simple set of parameters from which the computational grid is automatically generated. 1D and 2D grids for P2D and P3D models can also be generated.

![A selection of the parameterized battery geometries available. Clockwise from top left are: A single-layer pouch cell, CR 2016 coin cell, 30-layer pouch cell and jelly roll cylindrical cell. \label{fig:geometries}](figs/batterygeometries.pdf){width=100%}

# Graph based model development

Models grow in complexity as we extend or couple them. We introduce in BattMo a computational graph model design whose
aim is to help a developper to better control this complexity. A model is essentially defined through a computational
graph. The nodes are the variable names, as they can be introduced by the model, and should be recognized as such by a
user. The directed edges represents the functional relationships between the variables. A model is ready for simulation
when the root of the graph consists of the governing variables and the leaves as the governing equations. Interactive
tools to explore the computational graph of a model are available. Model hierarchy is an essential part of the
framework. Coupling two models is achieved by introducing a new coupling model, where the two models are sub-models. We
include their graph and add new edges corresponding to the coupling parts. We can focus on those as most part of the
sub-model remain unchanged. For more about this approach, see
[documentation](https://battmo.org/BattMo/computationalGraph/graphdoc.html ).

# Examples

Numerous documented and tested examples are provided with the code and demonstrates the features listed above. Multiple different data sets with different chemistries are also provided, including NMC, NCA and LNMO as well as Silicon enriched negative electrodes, standard electrolyte models, separators and current collector material properties. We seek to be cross-compatible, allowing for any chemistry and material combination with any battery format. A selection of the examples are documented at [https://battmo.org/BattMo/](https://battmo.org/BattMo/). The complete list of examples are available at [https://github.com/BattMoTeam/BattMo/tree/main/Examples](https://github.com/BattMoTeam/BattMo/tree/main/Examples).

# BattMo family

The following software include the BattMo family:

| Software                                             | Description                                                         |
|------------------------------------------------------|---------------------------------------------------------------------|
| [BattMo](https://github.com/BattMoTeam/BattMo)       | MATLAB and original version of BattMo presented in this publication |
| [BattMo.jl](https://github.com/BattMoTeam/BattMo.jl) | Julia version                                                       |
| [PyBattMo](https://github.com/BattMoTeam/PyBattMo)   | Python wrapper around BattMo.jl                                     |
| [BattMoApp](https://app.battmo.org/)                 | Web-application built on top of BattMo.jl                           |


# Acknowledgements

We acknowledge contributions from the European Union, Grant agreements 101069765 (IntelLiGent), 875527 (HYDRA), 957189 (BIG-MAP), 101104013 (BATMAX) and 101103997 (DigiBatt).

# References
