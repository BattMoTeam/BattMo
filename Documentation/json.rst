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

   The schemas may have common property fields. In the validation process, each schema is handled parallelly.

For example a component (say the electolyte) has its geometrical and material properties specified in two separate
schemas. We think it clarifies the presentation to separate the presentation of those properties and it corresponds also
to a convenient way to organize the input. We can easily switch between different geometrical models while keeping the
same material properties. For that, we use the :battmo:`mergeJsonStructs` function, see example
:ref:`here<mergeJsonStructs>`.

This is the schema for the input json file. It uses other separate schemas which takes care of specific parameter sets.

.. literalinclude:: ../Utilities/JsonSchemas/Simulation.schema.json
   :language: json

The schemas are available in the directory :battmofile:`JsonSchemas<Utilities/JsonSchemas>`. Here, we present the most
important ones.
              

Material Parameters
===================

The main schema is given by

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

.. literalinclude:: ../Utilities/JsonSchemas/Battery.schema.json
   :language: json

Electrode
---------

.. literalinclude:: ../Utilities/JsonSchemas/Electrode.schema.json
   :language: json

Coating
-------

.. literalinclude:: ../Utilities/JsonSchemas/Coating.schema.json
   :language: json


Interface
---------

.. literalinclude:: ../Utilities/JsonSchemas/Interface.schema.json
   :language: json


Solid Diffusion
---------------

.. literalinclude:: ../Utilities/JsonSchemas/SolidDiffusionModel.schema.json
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
                                                                    
