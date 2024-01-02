=========================
BattMo Model Architecture
=========================

We use a multi-model approach. The models are organized in a hierarchy, meaning that a given model can have
sub-models. A model corresponds to a given physical system and defines the functions and variables that will be needed
to assemble the discretized governing equations of the system.

For the simulation of Lithium-ion battery, we have at the top of the model hierarchy the the battery model (see
:ref:`schema <json:Material Parameters>` for the description of the standard input parameters). The battery model
contains the following sub-models:

* Negative electrode model :ref:`(schema) <json:Electrode>`
* Positive electrode model :ref:`(schema) <json:Electrode>`
* Electrolyte model :ref:`(schema) <json:Electrolyte>`
* Separator model :ref:`(schema) <json:Separator>`
* ThermalModel model :ref:`(schema) <json:Thermal Model>`
* Control model :ref:`(schema) <json:Simulation Control Parameters>`

.. figure:: img/cutbatterygraph.png
   :target: _images/cutbatterygraph.png

The **negative and positive electrodes** are instances of the same electrode model. The standard input parameters for
electrode model are given in its :ref:`schema <json:Electrode>`. The electrode model has two sub-models:

* Coating model :ref:`(schema) <json:Coating>`
* Current Collector model :ref:`(schema) <json:Current Collector>`


.. figure:: img/electrodegraph.png
   :target: _images/electrodegraph.png
   :width: 50%
   :align: center


The current collector model is optional. In particular, for the 1D model it is common to not include it.
           
The standard input parameters of the **coating model** are given in the associated :ref:`schema <json:Coating>`. The
coating model has three sub-models:

* Active material model :ref:`(schema) <json:Active Material>`
* Binder :ref:`(schema) <json:Binder>`
* Conductive additive model :ref:`(schema) <json:Conducting Additive>`

.. figure:: img/coatinggraph.png
   :target: _images/coatinggraph.png
   :width: 70%
   :align: center
   :class: with-border

           
In the case of a composite material, the coating model will have a different structure with two active material models.

.. _ArchitectureActiveMaterial:
   
The input parameters for the **Active Material** are described in the associated :ref:`schema <json:Active Material>`. The active material is
organized in two sub-models

* Interface :ref:`(schema) <json:Interface>`
* SolidDiffusion :ref:`(schema) <json:Solid Diffusion>`

.. figure:: img/activematerialgraph.png
   :target: _images/activematerialgraph.png
   :width: 50%
   :align: center

We have implemented two solid diffusion model, see :ref:`here <soliddiffusion:Solid Diffusion Models>`.

The **Control** :ref:`(schema) <json:Simulation Control Parameters>`, **Separator** :ref:`(schema) <json:Separator>` and
**Thermal** :ref:`(schema) <json:Thermal Model>` models do not have sub-models.
           











