
.. _runSEIActiveMaterial:

Particle simulation with SEI layer growth
--------------------------------------------------------------------
*Generated from runSEIActiveMaterial.m*



.. code-block:: matlab

  % clear the workspace and close open figures
  clear
  close all


Import the required modules from MRST
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
load MRST modules

.. code-block:: matlab

  mrstModule add ad-core mrst-gui mpfa


Setup the properties of Li-ion battery materials and cell design
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

  jsonstruct = parseBattmoJson(fullfile('ParameterData','ParameterSets','Safari2009','anode_sei.json'));
  
  % Some shorthands used for the sub-models
  ne    = 'NegativeElectrode';
  pe    = 'PositiveElectrode';
  am    = 'ActiveMaterial';
  sd    = 'SolidDiffusion';
  itf   = 'Interface';
  sei   = 'SolidElectrodeInterface';
  sr    = 'SideReaction';
  elyte = 'Electrolyte';
  
  inputparams = SEIActiveMaterialInputParams(jsonstruct);
  
  inputparams.(sd).N  = 10;
  inputparams.(sei).N = 10;


Setup the model
^^^^^^^^^^^^^^^

.. code-block:: matlab

  % We use a stand alone model for the particle
  inputparams.standAlone = true;
  
  % We initiate the model
  model = SEIActiveMaterial(inputparams);


Setup initial state
^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

  Nsd  = model.(sd).N;
  Nsei = model.(sei).N;
  
  % Initial concentration value at the electrode
  cElectrodeInit   = 0.75*model.(itf).saturationConcentration;
  % Initial value of the potential at the electrode
  phiElectrodeInit = 0;
  % Initial concentration value in the electrolyte
  cElectrolyte     = 5e-1*mol/litre;
  % Temperature
  T                = 298.15;
  
  % The following datas come from :cite:`Safari_2009`
  % Porosity of the SEI film
  epsiSEI     = 0.05;
  % Solvent concentration in the bulk of the electrolyte
  cECsolution = 4.541*mol/litre;
  % Solvent concentration in the SEI film
  cECexternal = epsiSEI*cECsolution;
  
  % We compute the OCP from the given data and use it to assign electrical potential in electrolyte
  initState.T = T;
  initState.(sd).cSurface = cElectrodeInit;
  initState = model.evalVarName(initState, {itf, 'OCP'});
  
  OCP = initState.(itf).OCP;
  phiElectrolyte = phiElectrodeInit - OCP;
  
  % From the values computed above we set the values of the initial state
  initState.E                = phiElectrodeInit;
  initState.I                = 0;
  initState.(sd).c           = cElectrodeInit*ones(Nsd, 1);
  initState.(sei).c          = cECexternal*ones(Nsei, 1);
  initState.(sei).cInterface = cECexternal;
  initState.(sei).delta      = 5*nano*meter;
  initState.R                = 0;
  
  % we set also static variable fields
  initState.T = T;
  initState.(itf).cElectrolyte   = cElectrolyte;
  initState.(itf).phiElectrolyte = phiElectrolyte;
  initState.(sr).phiElectrolyte  = phiElectrolyte;
  initState.(sei).cExternal      = cECexternal;


Setup schedule
^^^^^^^^^^^^^^

.. code-block:: matlab

  % Reference rate which roughly corresponds to 1 hour for the data of this example
  Iref = 1.3e-4*ampere/(1*centi*meter)^2;
  
  Imax = 1e1*Iref;
  
  total = 1*hour*(Iref/Imax);
  n     = 100;
  dt    = total/n;
  step  = struct('val', dt*ones(n, 1), 'control', ones(n, 1));
  
  % rampup value for the current function, see rampupSwitchControl
  tup = dt;
  srcfunc = @(time) rampupControl(time, tup, Imax);
  
  cmin = (model.(itf).guestStoichiometry0)*(model.(itf).saturationConcentration);
  control.stopFunction = @(model, state, state0_inner) (state.(sd).cSurface <= cmin);
  control.src = srcfunc;
  
  schedule = struct('control', control, 'step', step);


Setup non-linear solver
^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: matlab

  nls = NonLinearSolver();
  nls.errorOnFailure = false;
  
  model.nonlinearTolerance = 1e-5;


Run simulation
^^^^^^^^^^^^^^

.. code-block:: matlab

  model.verbose = true;
  [~, states, report] = simulateScheduleAD(initState, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);


Plotting
^^^^^^^^

.. code-block:: matlab

  set(0, 'defaulttextfontsize', 15);
  set(0, 'defaultaxesfontsize', 15);
  set(0, 'defaultlinelinewidth', 3);
  set(0, 'defaultfigureposition', [10, 10, 800, 400]);
  
  ind = cellfun(@(state) ~isempty(state), states);
  states = states(ind);
  
  time = cellfun(@(state) state.time, states);
  
  cSurface = cellfun(@(state) state.(sd).cSurface, states);
  figure
  plot(time/hour, cSurface/(1/litre));
  xlabel('time / h');
  ylabel('Surface concentration / mol/L');
  title('Surface concentration');
  
  E = cellfun(@(state) state.E, states);
  figure
  plot(time/hour, E);
  xlabel('time / h');
  ylabel('Potential / V');
  title('Potential');
  
  
  cmin = cellfun(@(state) min(state.(sd).c), states);
  cmax = cellfun(@(state) max(state.(sd).c), states);
  
  for istate = 1 : numel(states)
      states{istate} = model.evalVarName(states{istate}, {sd, 'cAverage'});
  end
  
  caver = cellfun(@(state) max(state.(sd).cAverage), states);
  
  figure
  hold on
  plot(time/hour, cmin /(mol/litre), 'displayname', 'cmin');
  plot(time/hour, cmax /(mol/litre), 'displayname', 'cmax');
  plot(time/hour, caver/(mol/litre), 'displayname', 'total concentration');
  title('Concentration in particle / mol/L')
  legend show
  
  delta = cellfun(@(state) state.(sei).delta, states);
  figure
  plot(time/hour, delta/(nano*meter));
  xlabel('time [hour]');
  ylabel('thickness / nm');
  title('SEI thickness')
  
  c = states{end}.(sd).c;
  r = linspace(0, model.(sd).particleRadius, model.(sd).N);
  
  figure
  plot(r, c/(mol/litre));
  xlabel('radius / m')
  ylabel('concentration / mol/L')
  title('Particle concentration profile (last time step)')
  
  r = states{end}.(sei).delta;
  r = linspace(0, r, model.(sei).N);
  c = states{end}.(sei).c;
  
  figure
  plot(r/(nano*meter), c/(mol/litre));
  xlabel('x / mm')
  ylabel('concentration / mol/L');
  title('Concentration profile in SEI layer (last time step)');

.. figure:: runSEIActiveMaterial_01.png
  :figwidth: 100%

.. figure:: runSEIActiveMaterial_02.png
  :figwidth: 100%

.. figure:: runSEIActiveMaterial_03.png
  :figwidth: 100%

.. figure:: runSEIActiveMaterial_04.png
  :figwidth: 100%

.. figure:: runSEIActiveMaterial_05.png
  :figwidth: 100%

.. figure:: runSEIActiveMaterial_06.png
  :figwidth: 100%



complete source code can be found :ref:`here<runSEIActiveMaterial_source>`
