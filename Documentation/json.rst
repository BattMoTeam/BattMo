========================
JSON input specification
========================

The specification of the json input is provided through `JSON schemas
<https://json-schema.org/understanding-json-schema>`_. In the Json schema, we added extra properties that are not part
of the syntax but bring additional information such as a :code:`description` field and references to the onthology (work
in progress :todo:`add more here`).

The main Simulation schema contains schemas for

* A schema for the **geometry**
* A schema for the **battery cell material parameters**
* A schema for the **initialization of the state of the battery**
* A schema for the **time stepping**
* A schema for the **solver parameters**
* A schema for the **output specification**

.. note::

   The different schemas may have common object. In the validation process, each schema is handled parallelly. The
   function :battmofile:`mergeJsonStructs` can be used to compose different json files.

For example a component (say the electolyte) has its geometrical and material properties specified in two separate
schemas. Having separate schemas clarify the presentation and corresponds also to a convenient way to organize the
input. We can easily switch between different geometrical models while keeping the same material properties. For that,
we use the :battmo:`mergeJsonStructs` function, see the example :ref:`here<mergeJsonStructs>`.


Simulation Schema
=================

This is the main schema for the input json file. It uses other separate schemas which takes care of specific parameter sets.

.. literalinclude:: ../Utilities/JsonSchemas/Simulation.schema.json
   :language: json

The schemas are available in the directory :battmofile:`JsonSchemas<Utilities/JsonSchemas>`. Here, we present the most
important schemas.
              

Material Parameters
===================

The main schema for the material parameter is reproduced below
(:battmofile:`source<Utilities/JsonSchemas/Battery.schema.json>`). We recognize the battery model structure presented
:ref:`here<architecture:BattMo Model Architecture>`.

.. literalinclude:: ../Utilities/JsonSchemas/Battery.schema.json
   :language: json

It contains references to schemas that are written in separate files

* Electrolyte
* Electrode
  
  * Coating
    
    * Interface
    * Solid Diffusion
      
  * Current Collector
    
* Separator
* Thermal Model

We give the listing of those here.
              
Electrolyte
-----------

The ionic conductivity and the diffusion coefficient can be given as a function or a constant. When a function is given,
the json file should contain the function name that is used. The *signature* of the function is given in the schema.

.. literalinclude:: ../Utilities/JsonSchemas/Electrolyte.schema.json
   :language: json

Electrode
---------

.. literalinclude:: ../Utilities/JsonSchemas/Electrode.schema.json
   :language: json

Coating
-------

.. literalinclude:: ../Utilities/JsonSchemas/Coating.schema.json
   :language: json

Active Material
---------------

.. literalinclude::  ../Utilities/JsonSchemas/ActiveMaterial.schema.json
   :language: json
              
Interface
---------

.. literalinclude:: ../Utilities/JsonSchemas/Interface.schema.json
   :language: json

Solid Diffusion
---------------

.. literalinclude:: ../Utilities/JsonSchemas/SolidDiffusionModel.schema.json
   :language: json

Binder
------

.. literalinclude::  ../Utilities/JsonSchemas/Binder.schema.json
   :language: json
              
Conducting Additive
-------------------

.. literalinclude::  ../Utilities/JsonSchemas/ConductingAdditive.schema.json
   :language: json

Current Collector
-----------------

.. literalinclude:: ../Utilities/JsonSchemas/CurrentCollector.schema.json
   :language: json



Separator
---------

.. literalinclude:: ../Utilities/JsonSchemas/Separator.schema.json
   :language: json


Thermal Model
-------------

.. literalinclude:: ../Utilities/JsonSchemas/ThermalComponent.schema.json
   :language: json

                 
.. _geometryschema:                 

Geometry Setup
==============

.. literalinclude:: ../Utilities/JsonSchemas/Geometry.schema.json
   :language: json
                 
Simulation Control Parameters
=============================

.. literalinclude:: ../Utilities/JsonSchemas/ControlModel.schema.json
   :language: json


Time Stepping Parameters
========================

.. literalinclude:: ../Utilities/JsonSchemas/TimeStepping.schema.json
   :language: json


Solver Parameters
=================

.. literalinclude:: ../Utilities/JsonSchemas/Solver.schema.json
   :language: json

Output Parameters
=================

.. literalinclude:: ../Utilities/JsonSchemas/Output.schema.json
   :language: json
                                                                    
