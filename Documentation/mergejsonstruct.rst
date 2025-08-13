==================
Merging parameters
==================

Setup a P4D Model using :code:`mergeJsonStructs`
================================================

Let's build a pseudo-four-dimensional Li-ion battery model and explore more about how to mix-and-match |battmo| parameter definitions!

P4D simulations provide the most spatial detail for battery simulations. They can consider the effects of complex shapes like current collector tabs, jelly rolls, or multi-layer cells on overall performance. However, they are also the most computationally expensive simulations to run. They are therefore usually reserved for specicial cases when geometric effects play a significant role. 

P4D models can be setup using a similar JSON parameter definition as described in the section on Basic Usage.

|battmo| aims to provide a modular solution for building electrochemical models. That allows us to mix-and-match different sets of parameter definitions without needing to do much re-coding. In this example, we will use the :battmo:`mergeJsonStructs` function to build a model definition from multiple sources. 

Define Parameters
-----------------

We will combine five separate JSON files that define the parameters for:

- cell materials
- cell geometry
- control policy
- simulation settings
- output specifications

First, let's define our cell materials. We have provided a JSON file that contains material properties for an NMC-Graphite Li-ion battery and can parse this as a |battmo| structure:

.. code:: matlab
          
   jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', 'lithium_ion_battery_nmc_graphite.json');
   jsonstruct_material = parseBattmoJson(jsonfilename);

Next, we have defined the cell geometry properties in a separate JSON file that we can also parse into |battmo|:

.. code:: matlab
          
   jsonfilename = fullfile('Examples', 'JsonDataFiles', 'geometry3d.json');
   jsonstruct_geometry = parseBattmoJson(jsonfilename);            

We can take the same approach for the remaining parameters, as shown below, for:

Control policy

.. code:: matlab
          
   jsonfilename = fullfile('Examples', 'JsonDataFiles', 'cc_discharge_control.json');
   jsonstruct_control = parseBattmoJson(jsonfilename);         

Simulation settings

.. code:: matlab
          
   jsonfilename = fullfile('Examples', 'JsonDataFiles', 'simulation_parameters.json');
   jsonstruct_simparams = parseBattmoJson(jsonfilename);       

Output specifications

.. code:: matlab
          
   jsonfilename = fullfile('Examples', 'JsonDataFiles', 'extra_output.json');
   jsonstruct_output = parseBattmoJson(jsonfilename);         

We can now merge these parameter definitions into a single parameter set and run the simulation:

.. code:: matlab
          
   jsonstruct = mergeJsonStructs({jsonstruct_geometry , ...
                                  jsonstruct_material , ...
                                  jsonstruct_control  , ...
                                  jsonstruct_simparams, ...
                                  jsonstruct_output   , ...                               
                                 });

Run Simulation
--------------

.. code:: matlab
          
   output = runBatteryJson(jsonstruct);  

Visualize Results
-----------------

We plot the model using :battmo:`plotBatteryGrid` (note that the different axis are scaled differently)

.. code:: matlab
          
   model = output.model
   plotBatteryGrid(model)

.. figure:: img/3dmodel.png
   :target: _images/3dmodel.png
   
We find a extensive set of plotting functions in `MRST <https://www.sintef.no/Projectweb/MRST/>`_. You may be interested
to have a look at the `Visualization Tutorial
<https://www.sintef.no/projectweb/mrst/documentation/tutorials/visualization-tutorial/>`_. Let us use the
:mrstfile:`plotGrid<mrst-core/plotting/plotGrid.m>` and :mrstfile:`plotCellData<mrst-core/plotting/plotCellData.m>` to plot the
surface particle concentrations in both electrode at a given time step.
          
..
   The plots presented below are obtained using the script runExample3D_doc 

.. code:: matlab
          
   state = output.states{20};
   E = state.Control.E
   plotGrid(model.grid, 'facecolor', 'none', 'edgealpha', 0.1)
   plotCellData(model.NegativeElectrode.Coating.grid, state.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))
   plotCellData(model.PositiveElectrode.Coating.grid, state.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))
   title('Particle Surface Lithium Concentration');

.. figure:: img/3dconc.png
   :target: _images/3dconc.png
   

File links and insertions with :code:`parseBattmoJson`
======================================================

There are two mechanisms which can be used to combine JSON input files:

#. Merge function using :code:`mergeJsonStruct`
#. Direct insertion using :code:`parseBattmoJson`

We have just seen an example of the first mechanism, which can be used within Matlab when we setup the simulation.

The function :battmo:`parseBattmoJson` is used to parse a JSON input and create the corresponding matlab structure. It
basically relies on `jsondecode <https://se.mathworks.com/help/matlab/ref/jsondecode.html>`_.

In this process the reserved keyword properties :code:`isFile` combined with :code:`filename` are used to fetch and
insert in place JSON data located in separate files. Here is an example, taken from
:battmofile:`lithium_ion_battery_nmc_graphite.json<ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json>`
where we have the following lines

.. code:: json
          
  "NegativeElectrode": {
    "Coating": {
      "ActiveMaterial": {
        "Interface": {
          "isFile": true,
          "filename": "ParameterData/MaterialProperties/Graphite/graphite.json"
        }}}}

The content of the file :battmofile:`graphite.json<ParameterData/MaterialProperties/Graphite/graphite.json>` is then
inserted in place. Hence, when we write

.. code:: matlab

   filename = fileread('ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json')
   jsonstruct = parseBattmoJson(filename)

the :code:`jsonstruct` that is obtained is equivalent to the one where we would have copied and paste the content of
:battmofile:`graphite.json<ParameterData/MaterialProperties/Graphite/graphite.json>`.

.. collapse:: jsonstruct detail

   .. code:: json
             
     "NegativeElectrode": {
       "Coating": {
         "ActiveMaterial": {
           "Interface": {
             "saturationConcentration": 30555,
             "volumetricSurfaceArea": 723600,
             "density": 2240,
             "numberOfElectronsTransferred" : 1,
             "activationEnergyOfReaction": 5000,
             "reactionRateConstant": 5.031e-11,
             "guestStoichiometry100": 0.88551,
             "guestStoichiometry0": 0.1429,
             "chargeTransferCoefficient": 0.5,
             "openCircuitPotential" : {"type": "function",
             "functionname" : "computeOCP_Graphite_Torchio",
             "argumentlist" : ["cElectrode", "T", "cmax"]
             }}},          


