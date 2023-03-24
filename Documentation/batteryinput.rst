              


Battery cell
============

.. image:: img/cutbatterygraph.png
           :width: 100%
           :align: center

.. autoclass:: Battery.BatteryInputParams
   :members:
      

Electrode components
====================

An electrode (:class:`Electrode <Electrochemistry.Electrode>`) consists of an active material (:class:`ActiveMaterial
<Electrochemistry.ActiveMaterial>`), which contains an interface model (:class:`Interface <Electrochemistry.Interface>`)
and a solid diffusion model (:class:`SolidDiffusionModel <Electrochemistry.SolidDiffusionModel>`), and a current
collector (:class:`CurrentCollector <Electrochemistry.CurrentCollector>`)

.. image:: img/electrodegraph.png
           :width: 80%


Electrode
---------

.. autoclass:: Electrochemistry.ElectrodeInputParams
   :members:

      
Active Material
---------------

.. autoclass:: Electrochemistry.ActiveMaterialInputParams
   :members:

            
Interface
---------

.. autoclass:: Electrochemistry.InterfaceInputParams
   :members:

      
Solid Diffusion Models
----------------------

.. autoclass:: Electrochemistry.SolidDiffusionModelInputParams
   :members:

.. autoclass:: Electrochemistry.FullSolidDiffusionModelInputParams
   :members:

.. autoclass:: Electrochemistry.SimplifiedSolidDiffusionModelInputParams
   :members:

      
Current Collector
-----------------
      
.. autoclass:: Electrochemistry.CurrentCollectorInputParams
   :members:


Electrolyte components
======================

.. image:: img/electrolytegraph.png
           :width: 30%
           :align: center
                   
.. autoclass:: Electrochemistry.ElectrolyteInputParams
   :members:

.. autoclass:: Electrochemistry.SeparatorInputParams
   :members:
