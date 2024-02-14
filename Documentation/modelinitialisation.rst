============================
The Battery Simulation Model
============================

In the previous sections, we have seen how to run a simulation from a single json input. In this section, we detail the
simulation process. This is useful for a more advanced usage where a direct access to the solver is required, for
example for a case that is not covered by the json input interface.

To run a simulation, we need:

* A :battmo:`Battery` **model**, which knows how to setup and solve the governing equations of our problem,
* An **initial state** and
* A **schedule**, which provides the time stepping and can also contain setup for control (covered in :ref:`an other
  section<controlinput:Control model description>`).

We use the governing equations provided in the serie of papers :cite:`Latz2011,Latz2013,Latz2016`. They correspond to
the standard PXD model. The equations are **non-linear**. The charge and mass conservation equations are partial
differential equations. First, the equations are discretized. We use a finite volume method which ensures conservation
at the discrete level. For stability, we use a fully implicit method, also called backward Euler. For each time step, we
end us with a set of non-linear equations that we have to solve. We use Newton method.

The function :mrstfile:`simulateScheduleAD<mrst-autodiff/ad-core/simulators/simulateScheduleAD.m>` with the following
signature takes as argument an initial state, a model and a schedule and returns global variables, a cell array of
states and a report. Each state in the cell array corresponds to the solution computed at a given time step.

.. code:: matlab
          
   function [globVars, states, schedulereport] = simulateScheduleAD(initState, model, schedule, varargin)

In this function, the model has the task of assemblying the discrete residual equations, given the solution at the given
time step, and send them to a Newton solver.


Initialisation of a battery simulation model
============================================

To initialise a model using json input, we can use the function :battmo:`setupModelFromJson`. For example,

.. code:: matlab

   jsonstruct = parseBattmoJson('Examples/JsonDataFiles/sample_input.json')
   model = setupModelFromJson(jsonstruct)


Then, we obtain a :code:`model` which is an instance of :battmo:`Battery`.

Internally, the json input is converted into a **matlab input object** which we typically call :code:`inputparams` and
which is here an instance of :battmo:`BatteryInputParams`. An advanced user does not have to use the json interface but
can work only with the matlab input structures. This approach with a double layer for input (first json then matlab) can
appear redundant but the advantage is for the developper to operate within the same Matlab environment: The computation
models (for example :battmo:`Battery` model) is initiated using a matlab object (for this example
:battmo:`BatteryInputParams`). In particular, it enables us to run **validation** methods on the input (doing the same
operations directly on json files would have been significantly more complicated to implement).

Let us now have a closer look to the :code:`setupModelFromJson` function, line by line. Given a matlab structure
:code:`jsonstruct` obtained as above by parsing a json file, we first convert the data that is specified with units to
SI units,

.. code:: matlab

   jsonstruct = resolveUnitInputJson(jsonstruct);

.. note::
   
   Within the simulator, all the quantities are used with SI units.

Then, we construct the matlab input object:

.. code:: matlab
          
   inputparams = BatteryInputParams(jsonstruct);


To every submodel (see :ref:`architecture:BattMo Model Architecture` for an overview of those), there corresponds a
matlab input parameter object, which is given the name of the submodel with the suffix *InputParams*. For example,
corresponding to the :battmo:`ActiveMaterial` model, we find the :battmo:`ActiveMaterialInputParams` input
model. Typically, the property of the object corresponds to the property of the corresponding json input.

We add the geometry using the function :battmo:`setupBatteryGridFromJson` which uses the json input

.. code:: matlab

   inputparams = setupBatteryGridFromJson(inputparams, jsonstruct);          

Now the object :code:`inputparams` contains also the grids for each of the submodels. In general, the grids are
generated using so-called grid generator, see the :ref:`dedicated section<geometryinput:Battery Geometries>`. When we
use a standard parameterized geometry, then the grid parameters can be pass in the json structure (an example is given
:battmofile:`here<Examples/Documentation/jsonfiles/4680-geometry.json>` for a :ref:`Jelly Roll geometry<jellyroll>`). The
function :battmo:`setupBatteryGridFromJson` then takes care of calling the appropriate grid generator with the
parameter, as given in the json input file.

The input parameter object can now be validated. This step is important. In the :ref:`architecture:BattMo Model
Architecture`, sub-models can use the same parameters. However, the submodels are instantiate in a parallel manner. The
validate method which is called recursively at each model level can ensure the consistency of the submodels. (For
example in the :battmo:`ActiveMaterialInputParams<ActiveMaterialInputParams#105>`, we make sure that the :ref:`Interface
<ArchitectureActiveMaterial>` model and the :ref:`Solid Diffusion <ArchitectureActiveMaterial>` model use the same
volumetric surface area). We can thus call

.. code:: matlab
          
   inputparams = inputparams.validateInputParams();

and make sure our data is consistent. In fact, this function is called in the setup of model so that we do not need to
run it separately. Finally, we use our input parameter object to instantiate our battery model

.. code:: matlab

   model = Battery(inputparams)


Inspection of the model
=======================

The model contains all the parameters of the battery. You can inspect simply using the command window. There are
properties that are used by the solver that will remain obscure for a standard user. Yet, most of the names are explicit
enough and match with the :ref:`json schema definition<json:JSON input specification>` so that their meaning will be
clear. For example,

.. code:: matlab

   >> model
   
   model = 
   
     Battery with properties:
   
                              con: [1x1 PhysicalConstants]
                NegativeElectrode: [1x1 Electrode]
                PositiveElectrode: [1x1 Electrode]
                      Electrolyte: [1x1 Electrolyte]
                        Separator: [1x1 Separator]
                     ThermalModel: [1x1 ThermalComponent]
                          Control: [1x1 CCDischargeControlModel]
                              SOC: 0.9900
                            initT: 298.1500
   ...
   
   
Here, we recognize the :ref:`battery model architecture<architecture:BattMo Model Architecture>`. Just as an example, we
can look at the properties of the active material in the negative electrode


.. code:: matlab

   >> model.NegativeElectrode.Coating.ActiveMaterial

   ans = 
   
     ActiveMaterial with properties:
   
                    Interface: [1x1 Interface]
               SolidDiffusion: [1x1 FullSolidDiffusionModel]
       electronicConductivity: 100
                      density: 2240
                 massFraction: 0.9400
          thermalConductivity: 1.0400
         specificHeatCapacity: 632

and we recognize the property names and values given in the :battmofile:`input json
file<Examples/JsonDataFiles/sample_input.json#13>` that is used in this example.

Some properties of the model are computed at initialisation. This is the case for example of the effective electronic
conducitivities. Therefore,

.. warning::

   In general, you should never change the properties of the model directly. You can do so if you know the model in
   details.

   The reason is that some parameters are used to compute other dependent parameters. This computation is done at the
   model setup.

  
To change a model parameter, you can either do it in your json input structure or, as described earlier, using the
matlab input parameter object (:code:`inputparams`). The effectrive electronic conductivity of the coating in the
negative electrode is

.. code:: matlab

   >> model.NegativeElectrode.Coating
   
   ans = 
   
     Coating with properties:
   
                        ActiveMaterial     : [1x1 ActiveMaterial]
                        Binder             : [1x1 Binder]
                        ConductingAdditive : [1x1 ConductingAdditive]

                        ...
                        
                electronicConductivity     : 100.3328
       effectiveElectronicConductivity     : 82.5961          

                        ...

The intrinsic electronic conductivity of the coating is computed from the electronic conductivity of its constituent
(active material, binder, conducting additive). Then the *effective* electronic conductivity, which is used in the
charge conservation equation, takes into account the coating volume fraction and the Bruggeman coefficient.

Let us change the conductivity of the active material from 100 to 120 (remember we always use SI inside the code so that
those values are in Siemens/cm^2). We can proceed as follows

.. code:: matlab
          
   jsonstruct = parseBattmoJson('JsonDataFiles/sample_input.json');
   [model, inputparams] = setupModelFromJson(jsonstruct);
   % Change the value the electronic conductivity
   inputparams.NegativeElectrode.Coating.ActiveMaterial.electronicConductivity = 120;
   model = Battery(inputparams);

Then, 

.. code:: matlab

   >> model.NegativeElectrode.Coating
   
   ans = 
   
     Coating with properties:
   
                electronicConductivity: 118.4873
       effectiveElectronicConductivity: 97.5413

Computing and inspecting some standard static properties of the model
=====================================================================

For a battery cell, utility functions are available to compute the standard properties listed below

* :battmo:`computeCellMass` computes the **mass** of the battery and its components
* :battmo:`computeCellCapacity` computes the **capacity** of the the electrodes
* :battmo:`computeCellEnergy` computes the **total energy** of the battery when discharged at equilibrium conditions.
  It means that the transport effects are totally neglicted and corresponds to the case of an infinitly small CRate.

To print to screen all these properties, we can use conveniently an instance of the :battmo:`CellSpecificationSummary`
as shown below.

.. code:: matlab

   jsonstruct = parseBattmoJson('JsonDataFiles/sample_input.json');
   model = setupModelFromJson(jsonstruct);

   css = CellSpecificationSummary(model);

Then, using the :code:`printSpecifications` method, we get

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
   
We can also mention here the utility function :battmo:`computeCellEnergyGivenCrate`, even it is not a *static* property.
The function computes the energy produced by a cell for a given CRate.

.. code:: matlab

   >> output = computeCellEnergyGivenCrate(model, 2);
   >> fprintf('Energy at Crate=2 : %g [Wh]', output.energy / (watt*hour));
   
   Energy at Crate = 2 : 0.0110781 [Wh]

