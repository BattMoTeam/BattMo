================
Input Parameters
================

We need to propagate inputs in each of submodels. To do that we use a hierarchy of input classes that is sent
recursively to the models. The input class is named accordingly to the model it corresponds to and its source is located
in the same folder.

Input parameters for the electrode components
=============================================


An electrode (:class:`Electrode <Electrochemistry.Electrodes.Electrode>`) consists of an electrode active component
(:class:`ElectrodeActiveComponent <Electrochemistry.ElectrodeActiveComponent>`), which contains the active meterial
(:class:`ActiveMaterialInputParams <Electrochemistry.Electrodes.ActiveMaterialInputParams>`), and a current collector
(:class:`CurrentCollector <Electrochemistry.CurrentCollector>`)

Electrode
---------

.. autoclass:: Electrochemistry.Electrodes.ElectrodeInputParams
   :members:

Electrode Active Component
--------------------------

.. autoclass:: Electrochemistry.ElectrodeActiveComponentInputParams
   :members:

Active Material
---------------
      
.. autoclass:: Electrochemistry.Electrodes.ActiveMaterialInputParams
   :members:

Current Collector
-----------------
      
.. autoclass:: Electrochemistry.CurrentCollectorInputParams
   :members:


Input parameters for the Electrolyte components
===============================================

.. autoclass:: Electrochemistry.ElectrolyteInputParams
   :members:
.. autoclass:: Electrochemistry.SeparatorInputParams
   :members:
   
Input parameters for the Battery
================================

.. autoclass:: Battery.BatteryInputParams
   :members:

Input parameters for the Generic models
=======================================

.. autoclass:: Electrochemistry.ElectronicComponentInputParams
   :members:
.. autoclass:: Electrochemistry.ElectroChemicalComponentInputParams
   :members:
      
