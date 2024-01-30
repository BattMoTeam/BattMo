        
========================
JSON input specification
========================

The specification of the json input is provided through `JSON schemas
<https://json-schema.org/understanding-json-schema>`_. In the Json schema, we added extra properties that are not part
of the syntax but bring additional information such as a :code:`description` field and references to the onthology (see
`BattInfo <https://emmo-repo.github.io/domain-electrochemistry/index.html>`_ ).

The main Simulation schema contains schemas for

* A schema for the **geometry**
* A schema for the **battery cell material parameters**
* A schema for the **initialization of the state of the battery**
* A schema for the **time stepping**
* A schema for the **solver parameters**
* A schema for the **output specification**

.. note::

   The different schemas may have common object. In the validation process, each schema is handled parallelly. The
   function :battmo:`mergeJsonStructs` can be used to compose different json files.

For example a component (say the electolyte) has its geometrical and material properties specified in two separate
schemas. Having separate schemas clarify the presentation and corresponds also to a convenient way to organize the
input. We can easily switch between different geometrical models while keeping the same material properties. For that,
we use the :battmo:`mergeJsonStructs` function, see the example :ref:`here<mergeJsonStructs>`.

Here, we give an :ref:`Example<jsonexample:Json input example>`


Simulation Schema
=================

This is the main schema for the input json file. It uses other separate schemas which takes care of specific parameter sets.

.. literalinclude:: ../Utilities/JsonSchemas/Simulation.schema.json
   :language: json

The schemas are available in the directory
:battmofile:`JsonSchemas<Utilities/JsonSchemas>`. :ref:`Here<jsonexample:Simulation>` is an example of the main input
file where the inputs are given in separate files using the :code:`isFile` key.
              

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

See json :ref:`input example<jsonexample:Battery>`.  
              
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

See json :ref:`input example<jsonexample:Electrolyte>`.

Electrode
---------

The electrode input data contains essentially the input data for the coating and the current collector. The current collector is optional.

.. literalinclude:: ../Utilities/JsonSchemas/Electrode.schema.json
   :language: json

See json :ref:`input example<jsonexample:Negative Electrode>` and also :ref:`here<jsonexample:Positive Electrode>`.
              
Coating
-------

In the coating input data, we find input data for the active material, the binder and the conducting additive. In the
case where we have a composite material (:code:`active_material_type` is set to :code:`composite`), then we have to
provide the data for the two active materials (:code:`ActiveMaterial1` and :code:`ActiveMaterial1`).

.. literalinclude:: ../Utilities/JsonSchemas/Coating.schema.json
   :language: json

See json :ref:`input example<jsonexample:Coating>`

Active Material
---------------

The active material input data contains the conductivity, the thermal data and the input data for the interface and
solid diffusion. The interface refers to the processes that occur there, i.e. the chemical reactions, see below. The
property :code:`diffusionModelType` is used to choose between the different diffusion model available, see also
:ref:`here<soliddiffusion:Solid Diffusion Models>`.

.. note::

   The *model switch* for the diffusion model (i.e. :code:`diffusionModelType`) is provided in the model *above* the
   diffusion model itself, in this case the active material model. When we initialise a sub-model, we need to know its
   type. By having the *model switch* in the model above, we can directly choose and start the corresponding
   initializaton. This design choice is in fact used consistently in BattMo.
   
.. literalinclude::  ../Utilities/JsonSchemas/ActiveMaterial.schema.json
   :language: json

See json :ref:`input example<jsonexample:Active Material>`

Interface
---------

The interface input data gives the specification of the chemical reaction occuring there. In particular, we find the
definition of open circuit potential (:code:`openCircuitPotential`). As mentioned
:ref:`above<note-on-future-function-support>`, we plan to include support for tabulated and string input for
functions.

.. literalinclude:: ../Utilities/JsonSchemas/Interface.schema.json
   :language: json
   
See json :ref:`input example<jsonexample:Interface>`
              
Solid Diffusion
---------------

The solid diffusion input data contains in particular the particle radius, the volumetric surface area and the reference
diffusion coefficient. We use an Arrhenius-type of equation to compute the diffusion coefficient. In the case of the
full diffusion model, we can provide a diffusion coefficient that depends on the concentration, see below.

.. literalinclude:: ../Utilities/JsonSchemas/SolidDiffusionModel.schema.json
   :language: json

See json :ref:`input example<jsonexample:Solid Diffusion>`
              
Full Solid Diffusion
--------------------

We can specify the diffusion coefficient as a function of the state of charge. To compute the state of charge, we need
the guest stoichiometries and the saturation concentration.

.. note::

   The guest stoichiometries and the saturation concentration are also input data for the
   :ref:`interface<json:Interface>`. The given values should then be the same. However, the submodels are initiated
   parallelly at the model level. For that reason, we use a validation mechanism that checks the consistency of the data
   of the submodels, when needed, see :battmo:`here <ActiveMaterialInputParams#107>` for an example.


.. literalinclude:: ../Utilities/JsonSchemas/FullSolidDiffusionModel.schema.json
   :language: json
              
See json :ref:`input example<jsonexample:Solid Diffusion>`
              
Binder
------

The conductivity of the binder and conductiving additive are used to compute the overall conductivity of the coating.

.. literalinclude::  ../Utilities/JsonSchemas/Binder.schema.json
   :language: json

See json :ref:`input example<jsonexample:Binder>`

Conducting Additive
-------------------

The conductivity of the binder and conductiving additive are used to compute the overall conductivity of the coating.

.. literalinclude::  ../Utilities/JsonSchemas/ConductingAdditive.schema.json
   :language: json

See json :ref:`input example<jsonexample:Conducting Additive>`
              
Current Collector
-----------------

.. literalinclude:: ../Utilities/JsonSchemas/CurrentCollector.schema.json
   :language: json

See json :ref:`input example<jsonexample:Current Collector>`

Separator
---------

The porosity of the separator gives the volume fraction of the electrolyte in that region.

.. literalinclude:: ../Utilities/JsonSchemas/Separator.schema.json
   :language: json

See json :ref:`input example<jsonexample:Separator>`


Thermal Model
-------------

The thermal parameters such as thermal capacity and conductivity are part of the material parameters. In the thermal
model, we include the external temperature and the heat transfer paremeters with the exterior domain. The later depend
often on the geometry, and they are in fact also included in the schema there, see below. We have included a flag to
indicate if we consider wet or dry properties. This flag is not yet supported and we always consider dry properties,
from which the effective wet properties are computed.

.. literalinclude:: ../Utilities/JsonSchemas/ThermalComponent.schema.json
   :language: json

See json :ref:`input example<jsonexample:Thermal Model>`

.. _geometryschema:                 

Geometry Setup
==============

BattMo supports very general grid structures. However, the geometrical models must be *constructed*. Typically, in the
manufacturing industry, this operation is done using `CAD software
<https://en.wikipedia.org/wiki/Computer-aided_design>`_. For batteries, there exist standard designs which can be
parameterized by a small set of parameters, see the :ref:`dedicated page<geometryinput:Battery Geometries>`.

For each design, the parameters are described in the schema.

.. literalinclude:: ../Utilities/JsonSchemas/Geometry.schema.json
   :language: json

See json :ref:`input example<jsonexample:Geometry>`

                 
Simulation Control Parameters
=============================

The control options are presented :ref:`here<controlinput:Control models>`. A description of each parameters for the
various control models can be read from the schema.

.. literalinclude:: ../Utilities/JsonSchemas/ControlModel.schema.json
   :language: json

See json :ref:`input example<jsonexample:Control>`
              

Time Stepping Parameters
========================

The description of the time stepping parameters can be read from the schema. Default parameters depending on the chose control model are provided.

.. literalinclude:: ../Utilities/JsonSchemas/TimeStepping.schema.json
   :language: json

See json :ref:`input example<jsonexample:Time Stepping>`
              
Solver Parameters
=================

Default parameters for the solver are provided. There exist a json interface to modify those and the corresponding
parameters are desribed in the schema. Many more options are available at the matlab level, which we do not document
here .

.. literalinclude:: ../Utilities/JsonSchemas/Solver.schema.json
   :language: json

Output Parameters
=================

Some extra post-processed parameters can be asked for already at the json level. Otherwise, it is always possible to
compute those afterwards (see function :battmo:`computeEnergyDensity` for example).

.. literalinclude:: ../Utilities/JsonSchemas/Output.schema.json
   :language: json
                                                                    
