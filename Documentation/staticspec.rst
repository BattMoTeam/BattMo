===================
Cell Specifications
===================

For a given battery model, we can compute the mass, the capacity and the nominal energy without running any
simulation. A :code:`CellSpecificationSummary` object can be instantiated from a model to print all this information

We setup a model using :code:`runBatteryJson` but with the option :code:`runSimulation` set as false to avoid running the simulation.

.. code:: matlab

   jsonstruct = parseBattmoJson('JsonDataFiles/sample_input.json');
   model = runBatteryJson(jsonstruct, 'runSimulation', false);

   css = CellSpecificationSummary(model);

Then, using the :code:`printSpecifications` method, we get an overview of the cell specifications

.. code:: matlab

   >> css.printSpecifications
   
                  Packing mass : 0 [kg]
                   Temperature : 24.85 [C]
                          Mass : 3.59526e-05 [kg]
                        Volume : 1.36e-05 [L]
                Total Capacity : 0.00301148 [Ah]
   Negative Electrode Capacity : 0.00310324 [Ah]
   Positive Electrode Capacity : 0.00301148 [Ah]
                     N/P ratio : 1.03047 [-]
                        Energy : 0.0115753 [Wh]
               Specific Energy : 321.958 [Wh/kg]
                Energy Density : 851.122 [Wh/L]
               Initial Voltage : 4.17686 [V]   
   
There exist separate functions to compute all this information separatly.
               

* :battmo:`computeCellMass` computes the **mass** of the battery and its components
* :battmo:`computeCellCapacity` computes the **capacity** of the the electrodes
* :battmo:`computeCellEnergy` computes the **total energy** of the battery when discharged at equilibrium conditions.
  It means that the transport effects are totally neglicted and corresponds to the case of an infinitly small CRate.

   
