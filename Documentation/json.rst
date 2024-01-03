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
the json file should contain the function name that is used. The *signature* of the function is given in the schema in
form of a *argument list*. For example, below, we can read that that the :code:`ionicConductivity` is a function of
concentration and temperature, see examples :battmofile:`here<ParameterData/MaterialProperties/OrganicLiPF6Solutions>`.

.. _note-on-future-function-support:

.. note::

   At the moment, only function implemented in matlab are supported. To add a function, the user must therefore create
   it in matlab, with the right signature, and make sure it is in the path. We plan to add very soon a pure json
   interface, where a function can be given either as string (which is evaluated to obtain the value) or as a
   table. When given as a table, the value of the function is interpolated from the data points.

.. literalinclude:: ../Utilities/JsonSchemas/Electrolyte.schema.json
   :language: json

Electrode
---------

The electrode input data contains essentially the input data for the coating and the current collector. The current collector is optional.

.. literalinclude:: ../Utilities/JsonSchemas/Electrode.schema.json
   :language: json

Coating
-------

In the coating input data, we find input data for the active material, the binder and the conducting additive. In the
case where we have a composite material (:code:`active_material_type` is set to :code:`composite`), then we have to
provide the data for the two active materials (:code:`ActiveMaterial1` and :code:`ActiveMaterial1`).

.. literalinclude:: ../Utilities/JsonSchemas/Coating.schema.json
   :language: json

Active Material
---------------

The active material input data contains the conductivity, the thermal data and the input data for the interface and
solid diffusion. The interface refers to the processes that occur there, i.e. the chemical reactions.

.. literalinclude::  ../Utilities/JsonSchemas/ActiveMaterial.schema.json
   :language: json
              
Interface
---------

The interface input data gives the specification of the chemical reaction occuring there. In particular, we find the
definition of open circuit potential (:code:`openCircuitPotential`). As mentioned
:ref:`above<note-on-future-function-support>`, we plan to include support for tabulated and string data input for
functions.

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
                                                                    
