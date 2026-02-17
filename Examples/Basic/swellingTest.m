% BattMo Tutorial
% This tutorial explains how to setup and run a simulation in BattMo

%%% Setting up the environment
% BattMo uses functionality from :mod:`MRST <MRSTBattMo>`. This functionality 
% is collected into modules where each module contains code for doing 
% specific things. To use this functionality we must add these modules to 
% the matlab path by running:


mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

clear all
close all

% mrstDebug(20);

%%% Specifying the physical model
% In this tutorial we will simulate a lithium-ion battery consisting of a 
% negative electrode, a positive electrode and an electrolyte. *BattMo* 
% comes with some pre-defined models which can be loaded from JSON files.
% Here we will load the basic lithium-ion model JSON file which comes with
% Battmo.

fname = fullfile('ParameterData'        , ...
                 'BatteryCellParameters', ...
                 'LithiumIonBatteryCell', ...
                 'lithium_ion_battery_nmc_silicon.json');

jsonstruct = parseBattmoJson(fname);

%%%
% The parseBattmoJson function parses the JSON input and creates a matlab
% structure containing the same fields as the JSON input. This structure 
% can be changed to setup the model in the way that we want. 

%%%
% In this instance we will exclude temperature effects by setting
% use_thermal to false.

jsonstruct.use_thermal = false;

%%%
% We will also not use current collectors in this example:

jsonstruct.include_current_collectors = false;

%%%
% The structure created in the jsonstruct follows the same hierarchy as the
% fields in the JSON input file. These can be referenced by name in the
% jsonstruct. To make life easier for ourselves we define some shorthand
% names for various parts of the structure.

ne      = 'NegativeElectrode';
pe      = 'PositiveElectrode';
co      = 'Coating';
elyte   = 'Electrolyte';
thermal = 'ThermalModel';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';
cc      = 'CurrentCollector';

%%%
% Now we can set the diffusion model type for the active material (am) in the
% positive (pe) and negative (ne) electrodes to 'full'.

jsonstruct.(ne).(co).(am).diffusionModelType = 'swelling';
jsonstruct.(pe).(co).(am).diffusionModelType = 'full';

%%%
% To see which other types of diffusion model are available one can view 
% :class:`ActiveMaterialInputParams <Electrochemistry.ActiveMaterialInputParams>.

%%%
% When running a simulation, *BattMo* requires that all model parameters
% are stored in an instance of :class:`BatteryInputParams <Battery.BatteryInputParams>`. 
% This class is used to initialize the simulation and is accessed by
% various parts of the simulator during the simulation. This class is
% instantiated using the jsonstruct we just created:

inputparams = BatteryInputParams(jsonstruct);

%%%
% It is also possible to update the properties of this inputparams in a
% similar way to updating the jsonstruct. Here we set the discretisation
% level for the diffusion model. Other input parameters for the full diffusion
% model can be found here:
% :class:`FullSolidDiffusionModelInputParams <Electrochemistry.FullSolidDiffusionModelInputParams>`.

inputparams.(ne).(co).(am).(sd).N = 5;
inputparams.(pe).(co).(am).(sd).N = 5;

%%% Setting up the geometry
% Here, we setup the 1D computational mesh that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGenerator1D. Classes for generating other geometries can
% be found in the BattMo/Battery/BatteryGeometry folder.

gen = BatteryGeneratorP2D();

gen.xlength(2) = 1e-6;
gen.xlength(4) = 100e-6;
gen.sepnx = 10; % discretization number for negative current collector (default = 10)
gen.nenx  = 10; % discretization number for negative active material (default = 10)
gen.penx  = 10; % discretization number for separator (default = 10)

%%%
% Now, we update the inputparams with the properties of the mesh. This function
% will update relevent parameters in the inputparams object and make sure we have
% all the required parameters for the model geometry chosen.

inputparams = gen.updateBatteryInputParams(inputparams);

%%% Initialising the battery model object
% The battery model is initialized by sending inputparams to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.

model = GenericBattery(inputparams);
cgit = model.cgit;
cgit.printRootVariables
% model = model.setupComputationalGraph();
% cgt = model.computationalGraph;
% return

%%%
% In BattMo a battery model is actually a collection of submodels: 
% Electrolyte, Negative Electrode, Positive Electrode, Thermal Model and Control
% Model. The battery class contains all of these submodels and various other 
% parameters necessary to run the simulation.

%%% Plotting the OCP curves
% We can inspect the model object to find out which parameters are being
% used. For instance the information we need to plot the OCP curves for the
% positive and negative electrodes can be found in the interface structure
% of each electrode.

%T = 298.15;
%elde = {ne,pe};
%
%figure
%hold on
%for i = 1:numel(elde)
%    po_itf = model.(elde{i}).(am).(itf);
%
%    theta100 = po_itf.theta100;
%    theta0   = po_itf.theta0;
%    cmax     = po_itf.cmax;
%
%    soc   = linspace(0, 1);
%    theta = soc*theta100 + (1 - soc)*theta0;
%    c     = theta.*cmax;
%    OCP = po_itf.computeOCPFunc(c, T, cmax);
%
%    plot(soc, OCP)
%end
%xlabel('SOC [-]')
%ylabel('OCV [V]')
%title('OCV for both electrodes');
%legend(elde)

%%% Controlling the simulation
% The control model specifies how the simulation is controlled. This can
% also be thought of as the boundary conditions of the simulation.

%%%
% In the first instance we use IEswitch control policy.
% We set the total time scaled by the CRate in the model.
% The CRate has been set by the json file. We can access it here:

timestep.numberOfTimeSteps = 100;
% timestep.totalTime = 100*hour;

step    = model.Control.setupScheduleStep(timestep);
control = model.Control.setupScheduleControl();

% This control is used to set up the schedule
schedule = struct('control', control, 'step', step);


% tup = 0.1;
% switch model.Control.controlPolicy
%   case 'IEswitch'
%       switch model.Control.initialControl
%       case 'discharging'
%         inputI = model.Control.Imax;
%         inputE = model.Control.lowerCutoffVoltage;
%       case {'charging', 'first-charge'}
%         inputI = -model.Control.Imax;
%         inputE = model.Control.upperCutoffVoltage;
%       otherwise
%         error('initCase not recognized')
%     end
%     srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
%                                                 inputI, ...
%                                                 inputE);
%     control.IEswitch = true;
%     control = struct('src', srcfunc, 'IEswitch', true);
%   case 'CC'
%     srcfunc = @(time) rampupControl(time, tup, model.Control.Imax);
%     control.CC = true;
%   otherwise
%     error('control policity not recognized');
% end

%%%
% We create a control structure containing the source function and
% specifying that we want to use IESwitch control:

% control.src = srcfunc;

%%%
% Finally we collect the control and step structures together in a schedule
% struct which is the schedule which the simulation will follow:

% schedule = struct('control', control, 'step', step); 


%%% Setting the initial state of the battery
% To run simulation we need to know the starting point which we will run it
% from, in terms of the value of the primary variables being modelled at
% the start of the simulation. 
% The initial state of the model is setup using model.setupInitialState()
% Here we take the state of charge (SOC) given in the input and calculate
% equilibrium concentration based on theta0, theta100 and cmax.


initstate = model.setupInitialState(jsonstruct);

%%% Running the simulation
% Once we have the initial state, the model and the schedule, we can call
% the simulateScheduleAD function which will actually run the simulation:

%model.verbose = true;

nls = NonLinearSolver;
nls.errorOnFailure = false;

[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 


%% Plotting the results
% To get the results we use the matlab cellfun function to extract the
% values Control.E, Control.I and time from each timestep (cell in the cell
% array) in states. We can then plot the vectors.

ind = cellfun(@(x) ~isempty(x), states);
states = states(ind);

E = cellfun(@(x) x.Control.E, states); 
I = cellfun(@(x) x.Control.I, states);

T = cellfun(@(x) x.time, states); 

%% Plot E as a function of the time

figure()
tiledlayout

nexttile
plot(T/hour, E)
xlabel('time [hours]')
ylabel('Cell Voltage [V]')

nexttile
plot(T/hour, I)
xlabel('time [hours]')
ylabel('Cell Current [A]')


%% Plot the overpotential of the negative electrode as a function of time

for istate = 1 : numel(states)
    states{istate} = model.addVariables(states{istate});
end

%%

figure
hold on

negativeElectrodeSize = model.(ne).grid.cells.num;
L = {};

for i = 1 : negativeElectrodeSize

    eta = cellfun(@(state) state.(ne).(co).(am).(itf).eta(i), states);
    plot(T/hour, eta);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time / hours')
ylabel('\eta / V')
title('Overpotential in Negative Electrode')
legend(L);

%% plot of average concentration

figure
hold on

L = {};

for i = 1 : negativeElectrodeSize

    cAver = cellfun(@(state) state.(ne).(co).(am).(sd).cAverage(i), states);
    plot(T/hour, cAver);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time / hours')
ylabel('c / mol/m^3')
title('Average particle concentration')
legend(L);

%% plot of Porosity

figure
hold on

negativeElectrodeSize = model.(ne).grid.cells.num;
L = {};

for i = 1 : negativeElectrodeSize

    poro = cellfun(@(state) 1 - state.(ne).(co).volumeFraction(i), states);
    plot(T/hour, poro);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time / hour')
ylabel('porosity / 1')
title('Negative Electrode Porosity')
legend(L);

%% Plot the radius evolution

for istate = 1 : numel(states)
    states{istate} = model.evalVarName(states{istate}, {ne, co, am, sd, 'radius'});
end
    
figure
hold on

negativeElectrodeSize = model.(ne).grid.cells.num;
L = {};

for i = 1 : negativeElectrodeSize

    radius = cellfun(@(state) state.(ne).(co).(am).(sd).radius(i), states);
    plot(T/hour, radius);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time [hours]')
ylabel('radius / m')
title('Silicon particle radius')
legend(L);

%% plot volume fraction
    
figure
hold on

negativeElectrodeSize = model.(ne).grid.cells.num;
L = {};

for i = 1 : negativeElectrodeSize

    vf = cellfun(@(state) state.(ne).(co).volumeFraction(i), states);
    plot(T/hour, vf);
    L{end + 1} = "x = " + int2str(i);
    
end

xlabel('time [hours]')
ylabel('volume fraction')
title('volume fraction')
legend(L);

%% Total Lithium content

figure

m = [];
vols = model.(ne).(co).grid.cells.volumes;

for istate = 1 : numel(states)
    cMaxTot = model.(ne).(co).maximumTotalConcentration;

    state = states{istate};
    x     = state.(ne).(co).(am).(sd).x;
    
    m(end + 1) = sum(cMaxTot.*x.*vols);
    
end

F = PhysicalConstants.F;
plot(T/hour, m*F/hour);
xlabel('time [hours]')
ylabel('amount / Ah');
title('Lithium content in negative electrode')

%{
  Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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












