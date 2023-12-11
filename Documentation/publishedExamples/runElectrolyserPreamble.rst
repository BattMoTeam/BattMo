The implementation of a AEM simulation code was motivated by the need for better simulation capabilities for H2
electrolysis and the possibility to validate the code using an implementation we had access to, which is also the one
used in :cite:`gerhardt2023open`. The mathematical model for AEM includes many processes and we refer to the papers
mentioned above for all the details. Here, we include a figure from :cite:`gerhardt2023open` which contains the
following

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
robust independently of the time step size. We use automatic differentiation to assemble the equations and compute the
derivative of the system which is needed in the Newton algorithm. These computation steps are in fact generic and we
rely on the BattMo infrastructure to run them.

.. figure:: /img/electrolyser3.png  
   :target: ../_images/electrolyser3.png            
   :align: center
