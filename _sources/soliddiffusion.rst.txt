======================
Solid Diffusion Models
======================

We use the Doyler, Fuller and Newman model, see :cite:`Doyle_1993`. For each given spatial location, the concentration
profile in the solid is given for a representative spherical particle located there. The diffusion in the solid
electrode is then modeled as a diffusion equation in thee radial dimension of the pseudo particle. We use a standard
control volume method to solve the diffusion equation after discretization in the radial direction.

Reduced models for the solid diffusion have been considered to reduced the computation time. We have implemented the
*polynomial approximation method* as described in :cite:`Zhang_2007`.

The option :code:`diffusionModelType` in the Active Material part of the json schema (see :battmofile:`here
<Utilities/JsonSchemas/ActiveMaterial.schema.json#10>`) is used to choose the type of diffusion model.

In the script :battmo:`runChen2020`, we compare the two diffusion models in addition to the solution produced by `PyBaMM
<https://pybamm.org/>`_ for the Chen dataset :cite:`Chen2020DevelopmentModels`.

We reproduce the last plot produced in this script here. We observe a good overall match between the full and reduced
model for solid diffusion. The discripancies are large only at the beginning of the discharge. We observe also a good
agreement between BattMo and PyBaMM.

.. figure:: img/chencomp.png
   :target: _images/chencomp.png
   :width: 100%
   :align: center

   Discharge curves for Chen's dataset with the full and reduced model for diffusion. 


.. figure:: img/chencompzoom.png
   :target: _images/chencompzoom.png
   :width: 100%
   :align: center

   Zoom on the initial phase of the discharge. 



