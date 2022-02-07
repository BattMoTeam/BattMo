=============
Input classes
=============

We need to propagate inputs in each of submodels. To do that we use a hierarchy of input classes that is sent
recursively to the models. The input class is named accordingly to the model it corresponds to and its source is located
in the same folder.

Generic models
==============

.. autoclass:: Electrochemistry.ElectronicComponentInputParams
.. autoclass:: Electrochemistry.ElectroChemicalComponentInputParams

Electrode components
==========================================

.. autoclass:: Electrochemistry.ElectrodeActiveComponentInputParams
.. autoclass:: Electrochemistry.Electrodes.ActiveMaterialInputParams
.. autoclass:: Electrochemistry.CurrentCollectorInputParams
.. autoclass:: Electrochemistry.Electrodes.ElectrodeInputParams

Electrolyte components
======================

.. autoclass:: Electrochemistry.ElectrolyteInputParams
.. autoclass:: Electrochemistry.SeparatorInputParams
   
Battery
=======

.. autoclass:: Battery.BatteryInputParams

                        
