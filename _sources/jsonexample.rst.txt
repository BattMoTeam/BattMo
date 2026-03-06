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

Coating
=======

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/coating.json
   :language: json

Active Material
===============

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/activeMaterial.json
   :language: json

Interface
=========

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/interface.json
   :language: json
                            
Solid Diffusion
===============

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/solidDiffusion.json
   :language: json
                            
Binder
======

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/binder.json
   :language: json
                            
Conducting Additive
===================

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/conductingAdditive.json
   :language: json

Current Collector
=================

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/currentCollector.json
   :language: json

Positive Electrode
==================

Here, we provide all the parameter in a single json structure without linking to other files.

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/positiveElectrode.json
   :language: json

Electrolyte
===========

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/electrolyte.json
   :language: json

Separator
=========

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/electrolyte.json
   :language: json
              
Thermal Model
=============

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/thermalComponent.json
   :language: json

Geometry
========

We consider the 3D demo model, see :ref:`here<3dgeometry>`

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/geometry.json
   :language: json
              
Model Specification
===================

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/modelSpecification.json
   :language: json
              
Control
=======

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/control.json
   :language: json              
              
Time Stepping
=============

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/timeStepping.json
   :language: json
              
State Initialization
====================

.. literalinclude:: ../Examples/Documentation/jsonfiles/Example/stateInitialization.json
   :language: json                            
                                             
