===========
Basic Usage
===========

In this section, we describe how to run |battmo| using json inputs. The specification of json inputs is described in :ref:`json:Json input specification`

Standard Run
============

First, we load and parse the json input file using the commands

.. code:: matlab

   jsonstruct = parseBattmoJson('input.json')

Then, we run it as follows
  
.. code:: matlab

   output = runBatteryJson(jsonstruct)
                      
The json file :code:`input.json` given :ref:`here <jsoninputfile>` provides the needed data for the simulation. The
section :ref:`json:Json input specification` gives a complete description of each of the field. We have used long and
explicit names for a good readability.

The :code:`output` structure returns among other thing the model and the states. 

.. code:: matlab

   model  : [1x1 Battery]
   states : [1x1 struct]
          
We can plot the results using :battmo:`plotDashboard`. Here, for example at time step equal to 10,

.. code:: matlab

   plotDashboard(model, states, 'step', 10)


.. figure:: img/firstresult.png
   :target: _images/firstresult.png
   :width: 100%
   :align: center

   Dashboard for the solution at a given timestep.
   
     
Modifying the Json input directly from Matlab
=============================================

We can modify directly the json input by editing in the file. The json file is converted in a standard Matlab structure
using the Matlab in-built function `jsondecode <https://se.mathworks.com/help/matlab/ref/jsondecode.html>`_. We can
therefore modify it directly in matlab. Here, we modify the CRate values,

.. code:: matlab

   CRates = [0.8, 1, 2];
   for i = 1 : numel(CRates)
       jsonstruct.Control.CRate = CRates(i);
       output = runBatteryJson(jsonstruct);
       plotResult(output);
   end

For this example, we have writting a :code:`plotResult` function which extracts and plots from the output the time and
voltage values, see :ref:`here <plotResult>`.
   
.. figure:: img/crates.png
   :target: _images/crates.png
   :width: 70%
   :align: center   


Combining Json inputs
=====================


There are two mechanisms which can be used to combine json input files:

#. Direct insertion using :code:`parseBattmoJson`
#. Merge function using :code:`mergeJsonStruct`

Direct insertion using :code:`parseBattmoJson`
----------------------------------------------

The function :battmo:`parseBattmoJson` parses the json input to create the corresponding matlab structure, basically
relying on `jsondecode <https://se.mathworks.com/help/matlab/ref/jsondecode.html>`_. In this process the reserved
keyword properties :code:`isFile` combined with :code:`filename` are used to fetch and insert in place json data located
in separate files. Here is an example, taken from :battmofile:`lithium_ion_battery_nmc_graphite.json<ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json>` where we have the following lines

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
             "functionname" : "computeOCP_graphite",
             "argumentlist" : ["cElectrode", "T", "cmax"]
             }}},          

.. _mergeJsonStructs:

Merge function using :code:`mergeJsonStructs`
---------------------------------------------

We have implemented in Matlab a simple function that merge json files (feel free to implement it in your favorite
languages). The function :battmo:`mergeJsonStructs` takes a cell array of json structure parsed with `jsondecode
<https://se.mathworks.com/help/matlab/ref/jsondecode.html>`_ or :battmo:`parseBattmoJson` and merge the fields.

Let us look at an example where we change the geometry.  In :ref:`geometryinput:Battery Geometries`, we give an overview
of the various geometrical model we support.

We use the same material parameters as in the previous case,

.. code:: matlab
          
   jsonfilename = 'ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json';
   jsonstruct_material = parseBattmoJson(jsonfilename);

Let us consider the :code:`3d-demo` :ref:`case<3dgeometry>`. The 3D model can be found in the :battmofile:`Geometry
Schema<Utilities/JsonSchemas/Geometry.schema.json#113>`. We use the parameters given in
:battmofile:`geometry3d.json<Examples/jsondatafiles/geometry3d.json>` and fetch those using

.. code:: matlab
          
   jsonfilename = 'Examples/jsondatafiles/geometry3d.json';
   jsonstruct_geometry = parseBattmoJson(jsonfilename);            


We merge the two json inputs by calling

.. code:: matlab

   jsonstruct = mergeJsonStructs({jsonstruct_geometry , jsonstruct_material})


Now we have a json structure :code:`jsonstruct` that contains material properties and geometry obtained from two
separate files. After adding the rest of the simulation inputs as done in :battmo:`runJsonScript`, the simulation can be
run as before by running

.. code:: matlab

   output = runBatteryJson(jsonstruct);

We plot the model using :battmo:`plotBatteryMesh` (note that the different axis are scaled differently)

.. code:: matlab
          
   model = output.model
   plotBatteryMesh(model)

.. figure:: img/3dmodel.png
   :target: _images/3dmodel.png
   
We find a extensive set of plotting functions in `MRST <https://www.sintef.no/Projectweb/MRST/>`_. You may be interested
to have a look at the `Visualization Tutorial
<https://www.sintef.no/projectweb/mrst/documentation/tutorials/visualization-tutorial/>`_. Let us use the
:mrstfile:`plotGrid<core/plotting/plotGrid.m>` and :mrstfile:`plotCellData<core/plotting/plotCellData.m>` to plot the
surface particle concentrations in both electrode at a given time step.
          
..
   The plots presented below are obtained using the script runExample3D in Documentation/scripts/runExample3D

.. code:: matlab
          
   state = output.states{20};
   E = state.Control.E
   plotGrid(model.G, 'facecolor', 'none', 'edgealpha', 0.1)
   plotCellData(model.NegativeElectrode.Coating.G, state.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))
   plotCellData(model.PositiveElectrode.Coating.G, state.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface/(mol/litre))
   title('Particle Surface Lithium Concentration');

.. figure:: img/3dconc.png
   :target: _images/3dconc.png
   
