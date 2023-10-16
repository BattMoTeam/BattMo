
.. _runBattery1D:

Pseudo-Two-Dimensional (P2D) Lithium-Ion Battery Model
-------------------------------------------------------------------------
*Generated from runBattery1D.m*


This example demonstrates how to setup a P2D model of a Li-ion battery and run a simple simulation.

.. code-block:: matlab

  % Clear the workspace and close open figures
  clear all
  close all
  clc


Import the required modules from MRST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
load MRST modules

.. code-block:: matlab

  mrstModule add ad-core mrst-gui mpfa agmg linearsolvers


Setup the properties of Li-ion battery materials and cell design
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The properties and parameters of the battery cell, including the architecture and materials, are set using an instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`. This class is used to initialize the simulation and it propagates all the parameters throughout the submodels. The input parameters can be set manually or provided in json format. All the parameters for the model are stored in the paramobj object.

.. code-block:: matlab

  jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));
  
  % We define some shorthand names for simplicity.
  ne      = 'NegativeElectrode';
  pe      = 'PositiveElectrode';
  elyte   = 'Electrolyte';
  thermal = 'ThermalModel';
  am      = 'ActiveMaterial';
  itf     = 'Interface';
  sd      = 'SolidDiffusion';
  ctrl    = 'Control';
  cc      = 'CurrentCollector';
  
  jsonstruct.use_thermal = false;
  jsonstruct.include_current_collectors = false;
  
  jsonstruct.(pe).(am).diffusionModelType = 'full';
  jsonstruct.(ne).(am).diffusionModelType = 'full';
  
  jsonstruct.use_particle_diffusion = true;
  
  paramobj = BatteryInputParams(jsonstruct);
  
  paramobj.(ne).(am).InterDiffusionCoefficient = 0;
  paramobj.(pe).(am).InterDiffusionCoefficient = 0;
  
  paramobj.(ne).(am).(sd).N = 5;
  paramobj.(pe).(am).(sd).N = 5;
  
  paramobj = paramobj.validateInputParams();
  
  use_cccv = false;
  if use_cccv
      cccvstruct = struct( 'controlPolicy'     , 'CCCV',  ...
                           'initialControl'    , 'discharging', ...
                           'CRate'             , 1         , ...
                           'lowerCutoffVoltage', 2.4       , ...
                           'upperCutoffVoltage', 4.1       , ...
                           'dIdtLimit'         , 0.01      , ...
                           'dEdtLimit'         , 0.01);
      cccvparamobj = CcCvControlModelInputParams(cccvstruct);
      paramobj.Control = cccvparamobj;
  end


Setup the geometry and computational mesh
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Here, we setup the 1D computational mesh that will be used for the simulation. The required discretization parameters are already included in the class BatteryGenerator1D.

.. code-block:: matlab

  gen = BatteryGenerator1D();
  
  % Now, we update the paramobj with the properties of the mesh.
  paramobj = gen.updateBatteryInputParams(paramobj);


Initialize the battery model.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The battery model is initialized by sending paramobj to the Battery class constructor. see :class:`Battery <Battery.Battery>`.

.. code-block:: matlab

  model = Battery(paramobj);
  model.AutoDiffBackend= AutoDiffBackend();
  
  inspectgraph = false;
  if inspectgraph
      % plot the computational graph
      cgt = ComputationalGraphTool(model);
      cgt.getComputationalGraph('doplot', true);
      return
  end


Compute the nominal cell capacity and choose a C-Rate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The nominal capacity of the cell is calculated from the active materials. This value is then combined with the user-defined C-Rate to set the cell operational current.

.. code-block:: matlab

  CRate = model.Control.CRate;


Setup the time step schedule
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Smaller time steps are used to ramp up the current from zero to its operational value. Larger time steps are then used for the normal operation.

.. code-block:: matlab

  switch model.(ctrl).controlPolicy
    case 'CCCV'
      total = 3.5*hour/CRate;
    case 'IEswitch'
      total = 1.4*hour/CRate;
    otherwise
      error('control policy not recognized');
  end
  
  n  = 100;
  dt = total/n;
  step = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
  
  % we setup the control by assigning a source and stop function.
  % control = struct('CCCV', true);
  %  !!! Change this to an entry in the JSON with better variable names !!!
  
  switch model.Control.controlPolicy
    case 'IEswitch'
      tup = 0.1; % rampup value for the current function, see rampupSwitchControl
      srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                  model.Control.Imax, ...
                                                  model.Control.lowerCutoffVoltage);
      % we setup the control by assigning a source and stop function.
      control = struct('src', srcfunc, 'IEswitch', true);
    case 'CCCV'
      control = struct('CCCV', true);
    otherwise
      error('control policy not recognized');
  end
  
  % This control is used to set up the schedule
  schedule = struct('control', control, 'step', step);


Setup the initial state of the model
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The initial state of the model is setup using the model.setupInitialState() method.

.. code-block:: matlab

  initstate = model.setupInitialState();


Setup the properties of the nonlinear solver
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

  nls = NonLinearSolver();
  
  linearsolver = 'direct';
  switch linearsolver
    case 'agmg'
      mrstModule add agmg
      nls.LinearSolver = AGMGSolverAD('verbose', true, 'reduceToCell', false);
      nls.LinearSolver.tolerance = 1e-3;
      nls.LinearSolver.maxIterations = 30;
      nls.maxIterations = 10;
      nls.verbose = 10;
    case 'battery'
      nls.LinearSolver = LinearSolverBatteryExtra('verbose'     , false, ...
                                                  'reduceToCell', true, ...
                                                  'verbosity'   , 3    , ...
                                                  'reuse_setup' , false, ...
                                                  'method'      , 'direct');
      nls.LinearSolver.tolerance = 1e-4;
    case 'direct'
      disp('standard direct solver')
    otherwise
      error()
  end
  
  % Change default maximum iteration number in nonlinear solver
  nls.maxIterations = 10;
  % Change default behavior of nonlinear solver, in case of error
  nls.errorOnFailure = false;
  nls.timeStepSelector=StateChangeTimeStepSelector('TargetProps', {{'Control','E'}}, 'targetChangeAbs', 0.03);
  % Change default tolerance for nonlinear solver
  model.nonlinearTolerance = 1e-3*model.Control.Imax;
  % Set verbosity
  model.verbose = true;


.. code-block:: none

  standard direct solver



Run the simulation
^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

  [wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);


.. code-block:: none

  Solving timestep 001/100:                                            -> 50 Seconds, 399 Milliseconds
  =================================================================================================================================================================================================================================================================
  | It # | elyte_massCons (cell) | elyte_chargeCons (cell) | ne_am_chargeCons (cell) | pe_am_chargeCons (cell) | ne_am_sd_massCons (cell) | ne_am_sd_soliddiffeq (cell) | pe_am_sd_massCons (cell) | pe_am_sd_soliddiffeq (cell) | EIeq (cell) | controlEq (ctrl) |
  =================================================================================================================================================================================================================================================================
  |    1 |*0.00e+00              |*0.00e+00                |*0.00e+00                |*0.00e+00                |*0.00e+00                 |*0.00e+00                    |*0.00e+00                 |*0.00e+00                    | 3.11e+01    |*0.00e+00         |
  |    2 | 9.60e-01              | 5.73e-01                | 4.50e-01                | 2.25e-01                | 1.36e-01                 | 1.36e-01                    | 8.30e-02                 | 8.30e-02                    |*2.60e-09    |*0.00e+00         |
  |    3 |*2.77e-03              |*1.66e-03                |*1.36e-03                |*4.12e-04                |*4.09e-04                 |*4.09e-04                    |*1.52e-04                 |*1.52e-04                    |*2.60e-09    |*0.00e+00         |
  =================================================================================================================================================================================================================================================================
  ...


Process output and recover the output voltage and current from the output states.
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

  ind = cellfun(@(x) not(isempty(x)), states);
  states = states(ind);
  E = cellfun(@(x) x.Control.E, states);
  I = cellfun(@(x) x.Control.I, states);
  Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
  % [SOCN, SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
  time = cellfun(@(x) x.time, states);
  
  plot(time, E);

.. figure:: runBattery1D_01.png
  :figwidth: 100%


Plot the the output voltage and current
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
plotDashboard(model, states, 'step', 0);

.. code-block:: matlab

  %{


