=====================
BattMo Li-ion battery
=====================

We use a multi-model approach. The models are organized in a hierarchy, meaning that a given model can have
sub-models. A model is supposed to define all the functions and variables that will be needed to simulate a physical
system separately. At the top, we have a battery model with the submodels:

* Negative electrode model
* Positive electrode model
* Electrolyte model
* Separator model
* ThermalModel model
* Control model


.. figure:: img/cutbatterygraph.png
   :target: _images/cutbatterygraph.png


The negative and positive electrodes are instances of the same Electrode model. The electrode model has two sub-models:

* Coating model
* Current Collector model (optional)
  
The coating model has four sub-models:

* Active material model
* Binder
* Conductive additive model

In the case of a composite material, the coating model will have a different structure with two active material models (see )

.. toctree::
   :hidden: 

   Battery <batteryinput>












