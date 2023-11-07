========================
Json input specification
========================

The specification of the json input is provided through `JSON schemas
<https://json-schema.org/understanding-json-schema>`_. In the Json schema, we added extra properties that are not part
of the syntax but bring additional information such as a :code:`description` field and references to the onthology (work
in progress :todo:`add more here`).

The main Simulation schema contains schemas for

* A schema for the geometry
* A schema for the battery cell material parameters
* A schema for the initialization of the state of the battery
* A schema for the time stepping
* A schema for the solver parameters
* A schema for the output specification

.. note::

   The schemas may have common property fields. In the validation process, each schema is handled parallelly.

For example a component (say the electolyte) has its geometrical and material properties specified in two separate
schemas. We think it clarifies the presentation to separate the presentation of those properties and it corresponds also
to a convenient way to organize the input. We can easily switch between different geometrical models while keeping the
same material properties. For that, we use the :battmo:`mergeJsonStructs` function, see examples :todo:`add link`.

.. collapse:: Simulation Schema

   .. literalinclude:: ../Utilities/JsonSchemas/Simulation.schema.json
      :language: json


Geometry Setup
--------------

.. collapse:: Geometry Schema

   .. literalinclude:: ../Utilities/JsonSchemas/Geometry.schema.json
      :language: json


Material Parameters
-------------------

.. collapse:: Battery Schema

   .. literalinclude:: ../Utilities/JsonSchemas/Battery.schema.json
      :language: json

Simulation Control Parameters
-----------------------------

.. collapse:: Control Schema

   .. literalinclude:: ../Utilities/JsonSchemas/ControlModel.schema.json
      :language: json


Time Stepping Parameters
------------------------

.. collapse:: TimeStepping Schema

   .. literalinclude:: ../Utilities/JsonSchemas/TimeStepping.schema.json
      :language: json


Solver Parameters
-----------------

.. collapse:: Solver Schema
              
   .. literalinclude:: ../Utilities/JsonSchemas/Solver.schema.json
      :language: json

Output Parameters
-----------------

.. collapse:: Output Schema
              
   .. literalinclude:: ../Utilities/JsonSchemas/Output.schema.json
      :language: json
                                                                    
