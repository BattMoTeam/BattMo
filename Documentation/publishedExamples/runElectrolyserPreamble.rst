The implementation of a AEM simulation code was motivated by the need for better simulation capabilities for H2
electrolysis and the possibility to validate the code using an implementation we had access to, which is also the one
used in :cite:`gerhardt2023open`. The mathematical model for AEM is complex and includes many coupled processes . We
refer to the papers mentioned above for all the details. Here, we include a figure from :cite:`gerhardt2023open` which
illustrates the following processes

* Two parallel ionic conduction path ways: ionomer and liquid electrolyte,
* Water concentration gradient build-up in the membrane inducing water exchange between electrodes,
* Water exchange between the liquid and gas phases in the electrolyte,
* Hydroxide exchange between the electrolyte and ionomer.

.. figure:: /img/electolysermodeling.png  
   :target: ../_images/electolysermodeling.png
   :align: center
           
   Illustration taken from :cite:`gerhardt2023open`

We use a model hierarchy which enables us to structure the code in a hopefully clearly readable way. The main model is
split into three models

* Oxygen Evolution Electrode,
* Ionomer Membrane,
* Hydrogen Evolution Electrode.
  

.. figure:: /img/electrolyser1.png  
   :target: ../_images/electrolyser1.png
   :width: 100%
   :align: center

The electrode models share the same structure given by

* Porous Transport Layer,
* Catalyst Layer,
* Exchange Reactions Model,
           
.. figure:: /img/electrolyser2.png  
   :target: ../_images/electrolyser2.png
   :width: 90%
   :align: center
           
The governing equations consists of charge and mass conservations equations with exchange source terms, see Figure. They
are assembled at each model level for the core parts while the coupling terms are assembled at levels above in the model
hierarchy. We use generic discrete spatial differentiation operators that can be applied to any dimension, even if the
geometry we have setup for the result in this model is only 1D. Once the equations are discretized in space, we
discretize in time and solve the resulting non-linear system of equations using an implicit scheme. This approach is
robust independently of the choice oftime step size. We use automatic differentiation to assemble the equations and compute the
derivative of the system that are needed in the Newton algorithmused to solve
the equations. These computation steps are in fact generic and we
rely on the BattMo infrastructure to run them.

.. figure:: /img/electrolyser3.png  
   :target: ../_images/electrolyser3.png            
   :align: center

We consider a one-dimensional electrolyser cell where the volume fractions of each component are given in Figure below.
We increase the current linearly from 0 to 3 A/cm^2. The pressures at the boundary are set to 1 atmosphere and the
:math:`OH^{-}` concentration to 1 mol/litre.

.. figure:: /img/electrolyter_volume_fractions.png  
   :target: ../_images/electrolyter_volume_fractions.png
   :align: center

   Volume fractions used initially in the cell. Illustration picture taken from :cite:`gerhardt2023open`.

The input datas are given in a json format. A simple way to get acquainted to the format is to look at the examples used
in the example below. We provide also `json schemas <https://json-schema.org/>`_ that describes these inputs:

* Catalyst Layer  :battmofile:`CatalystLayer.schema.json <Electrolyser/JsonSchemas/CatalystLayer.schema.json>`
* Electrolyser  :battmofile:`Electrolyser.schema.json <Electrolyser/JsonSchemas/Electrolyser.schema.json>`
* Evolution Electrode  :battmofile:`EvolutionElectrode.schema.json <Electrolyser/JsonSchemas/EvolutionElectrode.schema.json>`
* Exchange Reaction  :battmofile:`ExchangeReaction.schema.json <Electrolyser/JsonSchemas/ExchangeReaction.schema.json>`
* Ionomer Membrane  :battmofile:`IonomerMembrane.schema.json <Electrolyser/JsonSchemas/IonomerMembrane.schema.json>`
* Porous Transport Layer  :battmofile:`PorousTransportLayer.schema.json <Electrolyser/JsonSchemas/PorousTransportLayer.schema.json>`
* Electrolyser Geometry  :battmofile:`ElectrolyserGeometry.schema.json <Electrolyser/JsonSchemas/ElectrolyserGeometry.schema.json>`
   
Here comes the listing of the code:


