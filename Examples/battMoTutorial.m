%% BattMo Tutorial
% This tutorial explains how to setup and run a simulation in BattMo

%%% Setting up the environment
% *BattMo* uses functionality from `MRST <MRSTBattMo>`. This functionality 
% is collected into modules where each module contains code for doing 
% specific things. To use this functionality we must add these modules to 
% the matlab path by running:

mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

%%% Specifying the physical model
% In this tutorial we will simulate a lithium-ion battery consisting of a 
% negative electrode, a positive electrode and an electrolyte. *BattMo* 
% comes with some pre-defined models which can be loaded from JSON files.
% Here we will load the basic lithium-ion model JSON file which comes with
% Battmo.
fname = fullfile('ParameterData','BatteryCellParameters',...
    'LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json');
jsonstruct = parseBattmoJson(fname);

% The parseBattmoJson function parses the JSON input and creates a matlab
% struct containing the same fields as the JSON input. This struct can be
% changed to setup the model in the way that we want. 

% In this instance we will exclude temperature effects by setting
% use_thermal to false.

jsonstruct.use_thermal = false;

% We will also not use current collectors in this example:

jsonstruct.include_current_collectors = false;

% Our model will simulate diffusion so we set use_particle_diffusion to
% true:

jsonstruct.use_particle_diffusion = true;

% The structure created in the jsonstruct follows the same hierarchy as the
% fields in the JSON input file. These can be referenced by name in the
% jsonstruct. To make life easier for ourselves we define some shorthand
% names for various parts of the structure.

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

% Now we can set the diffusion model type for the active material (am) in the
% positive (pe) and negative (ne) electrodes to 'full'.

jsonstruct.(pe).(am).diffusionModelType = 'full';
jsonstruct.(ne).(am).diffusionModelType = 'full';

% To see which other types of diffusion model are available one can view 
% <Electrochemistry.ActiveMaterialInputParams>.

% When running a simulation, *BattMo* requires that all model parameters
% are stored in an instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`. 
% This class is used to initialize the simulation and is accessed by
% various parts of the simulator during the simulation. This class is
% instantiated using the jsonstruct we just created:

paramobj = BatteryInputParams(jsonstruct);
paramobj = paramobj.validateInputParams();

% It is also possible to update the properties of this paramobj in a
% similar way to updating the jsonstruct. Here we set some more parameters
% for the diffusion model. The definitions for these are found in the 
% corresponding classes:
% :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>`
% and :class:`FullSolidDiffusionModelInputParams <Electrochemistry.FullSolidDiffusionModelInputParams>`.

paramobj.(ne).(am).InterDiffusionCoefficient = 0;
paramobj.(pe).(am).InterDiffusionCoefficient = 0;

paramobj.(ne).(am).(sd).N = 5;
paramobj.(pe).(am).(sd).N = 5;

%%% Setting up the geometry
% Here, we setup the 1D computational mesh that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGenerator1D. Classes for generating other geometries can
% be found in the BattMo/Battery/BatteryGeometry folder. 
gen = BatteryGenerator1D();

% Now, we update the paramobj with the properties of the mesh. This function
% will update relevent parameters in the paramobj object and make sure we have
% all the required parameters for the model geometry chosen.
paramobj = gen.updateBatteryInputParams(paramobj);

%%% Initialising the battery model object
% The battery model is initialized by sending paramobj to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
model = Battery(paramobj);

% In BattMo a battery model is actually a collection of submodels: 
% Electrolyte, Negative Electrode, Positive Electrode, Thermal Model and Control
% Model. The battery class contains all of these submodels and various other 
% parameters necessary to run the simulation.
% To see what properties the battery model object has we can print out the model 
% variable:
model

% We can also plot the computational graph using the ComputationalGraphTool in 
% BattMo.
cgt = ComputationalGraphTool(model);
cgt.getComputationalGraph('doplot', true);




%%% Controlling the simulation
% The control model specifies how the simulation is controlled. 
% In the first instance we use IEswitch control policy. (NOTE WHAT IS IESWITCH?)
% 












%%% Setting the initial state of the battery

% The initial state of the model is setup using the model.setupInitialState() method.

initstate = model.setupInitialState(); 




%%% Setting up solver properties


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



model.AutoDiffBackend= AutoDiffBackend();


%%% Running the simulation





[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 






%%% Plotting the results


ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
E = cellfun(@(x) x.Control.E, states); 
I = cellfun(@(x) x.Control.I, states);
Tmax = cellfun(@(x) max(x.ThermalModel.T), states);
% [SOCN, SOCP] =  cellfun(@(x) model.calculateSOC(x), states);
time = cellfun(@(x) x.time, states); 

plot(time, E);



%% Plot the the output voltage and current
% plotDashboard(model, states, 'step', 0);



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}












