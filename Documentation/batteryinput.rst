

Battery cell
============

A battery consists of two electrodes and an electrolyte. Each are implemented by a model which takes its own set of
input parameters. The convention for the input model class names is to write them using the model name followed by the
suffix :code:`InputParams`. In graph above, we can see that the battery model has submodel instances as properties. The property names are in black color while the model they belong to is in blue.

.. image:: img/cutbatterygraph.png
           :width: 100%
           :align: center

.. autoclass:: Battery.BatteryInputParams
   :members:
      

Electrode components
====================

An electrode consists of an active material, which contains an interface model and a solid diffusion model, and a
current collector. Each of those components have a own set of input parameters.

.. image:: img/electrodegraph.png
           :width: 80%


.. automodule:: Electrochemistry
                   
Electrode
---------

This model is derived from :class:`ElectronicComponentInputParams`

.. autoclass:: ElectrodeInputParams
   :members:

      
Active Material
---------------

.. autoclass:: ActiveMaterialInputParams
   :members:

            
Interface
---------

.. autoclass:: InterfaceInputParams
   :members:

      
Solid Diffusion Models
----------------------

.. autoclass:: SolidDiffusionModelInputParams
   :members:

.. autoclass:: FullSolidDiffusionModelInputParams
   :members:

.. autoclass:: SimplifiedSolidDiffusionModelInputParams
   :members:

      
Current Collector
-----------------
      
.. autoclass:: CurrentCollectorInputParams
   :members:


Electrolyte components
======================

.. image:: img/electrolytegraph.png
           :width: 30%
           :align: center
                   
.. autoclass:: ElectrolyteInputParams
   :members:

.. autoclass:: SeparatorInputParams
   :members:


Electronic Component
--------------------

Base model for all component with a electrical potential and a charge conservation equation

.. autoclass:: ElectronicComponentInputParams
   :members:
