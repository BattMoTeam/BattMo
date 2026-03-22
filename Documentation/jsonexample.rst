===================
My First Json Input
===================

In this page, we present example of json inputs. The parameter names haven hopefully chosen in a way that they will be
easily understood by battery modelers. Their description is included in the json schemas, see some links include below,
and this documentation :ref:`this page<json:JSON input specification>`.

The json input follows the battery model hierarchy as described in :ref:`BattMo Model Architecture<architecture:BattMo
Model Architecture>`. The input contains also simulation specific parameters, such as solver parameters and time
stepping.

The keywords :code:`isFile`, :code:`filename` and :code:`filenames` are recognized by the parser
(:code:`parseBattmoJson`) and used to include inputs from other files, as shown in the next example

Simulation
==========

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/simulation.json
   :language: json

The keyword :code:`isFile` is recognized. The value of :code:`filenames` contain an array with a list of the files that
we want to include. Those files are also the ones that are presented in this page (all the files are in
:battmofile:`this directory <Examples/Documentation/jsonfiles/Example/>`, which comes with your battmo installation)

The path of the file is first tested as such but if no file is found, the parser will check the path that is obtained
relative to your battmo installation path.

You can actually run this illustration example using

.. code:: matlab

   jsonstruct = parseBattmoJson('Examples/Documentation/jsonfiles/Example/simulation.json');
   output = runBatteryJson(jsonstruct);

Battery
=======

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/battery.json
   :language: json

Here, we recognize the battery hierarchy where the sumodels are the negative and positive electrodes, the electrolyte,
the separator and the thermal model, if thermal effects are included. The initial temperature and state of charge are
also included here.

Negative Electrode
==================

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/negativeElectrode.json
   :language: json

The electrode consists of a current collector and a coating. They are fetched from files. Here, we show only the
negative electrode.

Coating
=======

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/coating.json
   :language: json

The coating consists of the active material, the binder and the conducting additive. The coating has also effective
parameters that are set at this level, such as the :code:`effectiveDensity` and the :code:`bruggemanCoefficient`. The
effective density is used to compute the porosity from the properties of the three coating components.

Active Material
===============

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/activeMaterial.json
   :language: json

The submodels in the active material are the **interface** which collects the parameters related to the surface reaction
between the electode and the electrolyte and the **solid diffusion** model which collects the the parameters related to
the solid diffusion in the pseudo particle.

Interface
=========

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/interface.json
   :language: json

We do not explain in detail here all the parameters. We hope the names are chosen so that they will be easily
recognizable by battery modellers. Again, we refer the :ref:`json schema <json:Interface>` for a reference, see the
:code:`description` key.

The most important parameter of the electrochemical reaction is certainly the **open circuit potential**. The open
circuit potential is a function. In this example, the function is implemented as a matlab function and we call
:code:`named function` this type of format.

Functional parameters enter the DFN model at several places and we have implemented a specific input for those which is
explained in details in :ref:`Functional Parameters <functioninterface:Functional parameters>`
              
Solid Diffusion
===============

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/solidDiffusion.json
   :language: json

The particle radius and the diffusion coefficient governs the time scale of the diffusion. The volumetric surface
parameter is an important part of the volume averaging approach the DFN model uses, and it is not easy to estimate.

Binder
======

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/binder.json
   :language: json

The binder contains intrinsic parameters that are used to compute the effective coating material parameters for thermal
and electronic conductivity but no other process are modeled at this level. 
              
Conducting Additive
===================

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/conductingAdditive.json
   :language: json

As binder, conductive additive contains intrinsic parameters that are used to compute the effective coating material parameters for thermal
and electronic conductivity but no other process are modeled at this level.
              
Electrolyte
===========

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/electrolyte.json
   :language: json

This is an other example where we see functional parameters, the ionic conductivity and the diffusion coefficient, which
are functions of concentration and temperature.

Current Collector
=================

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/currentCollector.json
   :language: json

Thermal effects and electrical current is modeled at the current collector.
              

Separator
=========

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/electrolyte.json
   :language: json

The binder enters the model only through thermal effects.              
              
Thermal Model
=============

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/thermalComponent.json
   :language: json

The thermal properties, namely the thermal capacity and conductivity, are effective properties which are computed from
the intrinsic properties of each components, as we have seen above.

Here, we find the parameters that characterize the heat transfer with the exterior. They are effective parameters,
except for the external temperature.
              
Geometry
========

We consider the 3D demo model, see :ref:`here<3dgeometry>`

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/geometry.json
   :language: json

We consider geometries that are parameterized. It is easy to consider an arbitrary design as long as the user can
provide a mesh with tags for the component. We have however not setup a json interface for such usage.

We can see that the geometrical parameters are dispatched in each component where they belong.

Again, we refer to json schema for an explanation of the parameters, see :ref:`Geometry Setup <json:Geometry Setup>`.

Control
=======

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/control.json
   :language: json              

Here, we use a discharge scenario with a given discarge rate and lower cutoff voltage. An overview with explanation of the available controls is given in :ref:`Control models <controlinput:Control models>`.
              
Model Specification
===================

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/modelSpecification.json
   :language: json

Under the model specification key, we include global model setups. In this case, we include thermal effects and the
current collector are included. 

State Initialization
====================

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/stateInitialization.json
   :language: json                            

The simulator requires an initial state for the system. The default choice is to provide a state of charge (soc) and set
constant concentration and potential values equal to their respective equilibrium values. We provide also the initial
electrolyte concentration at the electrolyte level.
              
Time Stepping
=============

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/timeStepping.json
   :language: json

This information is used by the simulator to setup the time stepping. Default values are also provided, see :ref:`Time Stepping Parameters <jsonexample:Time Stepping Parameters>`.
