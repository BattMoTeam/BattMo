====================
Model initialisation
====================




How to initialise the model
---------------------------

Initial concentration and potential in the model can be calculated based on equilibrium values for the given input paramaters. This can be done by calling initialiseState on a BatteryModel object:

.. code:: matlab

    model.initialiseState()

Alternatively a state structure contatining concentration and potential for all the relevant submodels can be given.

Calculating intial concentration
--------------------------------

Concentration is calculated using the intial SOC of the battery.

Calculating initial potential
-----------------------------
