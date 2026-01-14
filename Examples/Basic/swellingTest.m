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

mrstDebug(20);

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

inputparams.(ne).(co).(am).(sd).N = 2;
inputparams.(pe).(co).(am).(sd).N = 2;

%%% Setting up the geometry
% Here, we setup the 1D computational mesh that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGenerator1D. Classes for generating other geometries can
% be found in the BattMo/Battery/BatteryGeometry folder.

gen = BatteryGeneratorP2D();
gen.sepnx = 2; % discretization number for negative current collector (default = 10)
gen.nenx  = 2; % discretization number for negative active material (default = 10)
gen.penx  = 2; % discretization number for separator (default = 10)

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

timestep.timeStepDuration = 100;

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

[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

return


%%%
% The outputs from the simulation are:
% - wellSols: which provides the current and voltage of the battery at each 
% timestep. (This naming convention is a hangover from MRST where we model
% reservoir injection via injection wells).
% - states: which contains the values of the primary variables in the model
% at each timestep.
% - reports: which contains technical information about the steps used in
% the numerical solvers.


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
subplot(2,2,1)
plot(T/hour, E)
xlabel('time [hours]')
ylabel('Cell Voltage [V]')

%% Plot I as a function of the time

subplot(2,2,2)
plot(T/hour, I)
xlabel('time [hours]')
ylabel('Cell Current [A]')

%% Plot the overpotential of the negative electrode as a function of time

for istate = 1 : numel(states)
    states{istate} = model.addVariables(states{istate});
end

%%

subplot(2,2,3)
hold on
negativeElectrodeSize = model.(ne).grid.cells.num;
L = "x = 1";

for i = 1 : negativeElectrodeSize

    eta = cellfun(@(state) state.(ne).(co).(am).(itf).eta(i), states);
    plot(T/hour, eta);

    if i > 1
        
        L(end+1) = "x = " + int2str(i);
        
    end
    
end
xlabel('time [hours]')
ylabel('Eta of the Negative electrode')
legend(L);

%% Plot of the porosity as a function of the time for different position across the negative electrode

% Plot the porosity as a function of the position for different times
subplot(2,2,4)
hold on
negativeElectrodeSize = gen.xlength(2);
N_elements_ne = gen.nenx;
deltaX = (negativeElectrodeSize/(N_elements_ne-1)) * 10^6;
totalTime = length(T);

position = [];
for x = 1:N_elements_ne
    position(end+1) = (x - 1)*deltaX;
end

legendTime = "t = 0 hour";
porosity = [];
for i = 1:N_elements_ne
    porosity(end + 1) = initstate.(ne).(co).porosity(i);
end

plot(position, porosity);

%% Plot of porosity as a function of time

for t = 1:totalTime
    hold on
    %Only draw the curve for timestep multiples
    if t < 80
        timestep = 8;
    else
        timestep = 110;
    end
    if mod(t,timestep) == 0
        porosity = [];
        for i = 1 : N_elements_ne
            porosity(end+1) = states{t}.(ne).(co).porosity(i);
        end
        
        plot(position, porosity);

        if t > 1
            t = T(t)/hour;
            legendTime(end+1) = "t = " + num2str(t,2) + " hour";
        end
    end 
end
xlabel('Position across the Negative Electrode (in µm)')
ylabel('Porosity of the Negative Electrode')
legend(legendTime);


%% Plot the porosity near the current collector on the negative electrode as a function of the state of charge

figure
 
ne_itf = model.(ne).(co).(am).(itf);
ne_sd  = model.(ne).(co).(am).SolidDiffusion;

totalTime = states{numel(states)}.time;

theta100 = ne_itf.guestStoichiometry100;
theta0   = ne_itf.guestStoichiometry0;
cmax     = ne_itf.saturationConcentration;
r0       = ne_sd.particleRadius;
F        = model.con.F;

soc      = [];
porosity = [];

for ind = 1 : numel(T)
    
    cAverage = states{ind}.(ne).(co).(am).(sd).cAverage(1);
    
    theta = cAverage/cmax;
    
    poros = states{ind}.(ne).(co).porosity(1);
    
    stoc = (theta - theta0)./(theta100 - theta0);
    
    soc(end + 1)      = stoc;
    porosity(end + 1) = poros;
    
end

plot(soc, porosity,'LineWidth',1.2);
axis([0 1 0 1]);
axis square;

xlabel('State of Charge near the current collector')
ylabel('Porosity near the current collector')

%% Plot the radius evolution

figure

ne_itf = model.(ne).(co).(am).(itf);
ne_sd  = model.(ne).(co).(am).SolidDiffusion;

Y = [];
X = [];


N = model.(ne).(co).(am).(sd).N;

for t = 1:numel(states)

    sumTheta = 0;
    for x = 1:N_elements_ne       
        sumConcentrations = 0;
        for i = 1:N
            c = states{t}.(ne).(co).(am).(sd).c((x-1)*N +i);
            sumConcentrations = sumConcentrations + c;
        end
        cAverage = sumConcentrations/N;
        
        radius = computeRadius(cAverage, cmax,r0);
    end
    
    Y(end+1) = radius;
    X(end+1) = T(t)/hour;
end

plot(X, Y);
ylabel('Silicon particle radius')
xlabel('Time (in hours)')


%% Plot SOC for each electrode as a function of time

figure
subplot(2,1,1)

ne_itf = model.(ne).(co).(am).(itf);
ne_sd  = model.(ne).(co).(am).SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

guestStoichiometry100 = ne_itf.guestStoichiometry100;
guestStoichiometry0   = ne_itf.guestStoichiometry0;
cmax                  = ne_itf.saturationConcentration;
r0                    = ne_sd.particleRadius;
F                     = model.con.F;
N                     = ne_sd.N;
N_elements_ne         = gen.nenx;


Y = [];
X = [];

for t = 1:numel(states)


    state = states{t};
    state.(ne).(am).SolidDiffusion = model.(ne).(co).(am).(sd).updateAverageConcentration(state.(ne).(co).(am).SolidDiffusion);
    cAverage = state.(ne).(co).(am).(sd).cAverage;

    vols = model.(ne).(co).grid.cells.volumes;

    theta = (sum(vols .* cAverage))./(sum(vols).* cmax);

    %sumTheta = 0;
    %    for x = 1:N_elements_ne       
    %    sumConcentrations = 0;
    %    for i = 1:N
    %        c = states{t}.(ne).(co).(am).(sd).c((x-1)*N +i);
    %        sumConcentrations = sumConcentrations + c;
    %    end
    %    cAverage = sumConcentrations/N;
    %
    %    theta = cAverage/cmax;
    %
    %    sumTheta = sumTheta + theta;
    %    end
    %    theta = sumTheta/N_elements_ne;
    %
    soc = (theta - guestStoichiometry0)/(guestStoichiometry100 - guestStoichiometry0);


    Y(end+1) = soc;
    X(end+1) = T(t)/hour;
end

plot(X, Y);

ylabel('State of Charge in the NEGATIVE ELECTRODE')
xlabel('Time (in hours)')


subplot(2,1,2)

po_itf = model.(pe).(co).(am).(itf);
po_sd  = model.(pe).(co).(am).SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

guestStoichiometry100 = po_itf.guestStoichiometry100;
guestStoichiometry0   = po_itf.guestStoichiometry0;
cmax                  = po_itf.saturationConcentration;
r0                    = po_sd.particleRadius;
F                     = model.con.F;
N                     = po_sd.N;
N_elements_pe         = gen.penx;


Y = [];
X = [];

for t = 1:totalTime

    sumTheta = 0;
    for x = 1:N_elements_pe
        sumConcentrations = 0;
        for i = 1:N
            c = states{t}.(pe).(co).(am).(sd).c((x-1)*N + i);
            sumConcentrations = sumConcentrations + c;
        end
        cAverage = sumConcentrations/N;    
        theta = cAverage/cmax;

        sumTheta = sumTheta + theta;
    end
    theta = sumTheta/N_elements_pe;

    soc = (theta-guestStoichiometry100)/(guestStoichiometry0-guestStoichiometry100);

    Y(end+1) = soc;
    X(end+1) = T(t)/hour;
end

plot(X, Y);

ylabel('State of charge in the POSITIVE ELECTRODE')
xlabel('Time (in hours)')

return

%%

figure

ne_itf = model.(ne).(co).(am).(itf);
ne_sd  = model.(ne).(co).(am).SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

guestStoichiometry100 = ne_itf.guestStoichiometry100;
guestStoichiometry0   = ne_itf.guestStoichiometry0;
cmax     = ne_itf.cmax;
r0       = ne_sd.rp;
F        = model.con.F;
N        = ne_sd.N;
N_elements_ne = gen.nenx;


Y = [];
X = [];

for t = 1:totalTime

    vf_ne = model.(ne).(co).(am).Interface.volumeFraction;
    vf_pe = model.(pe).(co).(am).Interface.volumeFraction;
    porosElyte = model.(elyte).volumeFraction;
    rp0_ne = model.(ne).(co).(am).(sd).rp;
    rp0_pe = model.(pe).(co).(am).(sd).rp;

    state = states{t};

    state.(ne).(co).(am).SolidDiffusion = model.(ne).(co).(am).(sd).updateAverageConcentration(state.(ne).(co).(am).SolidDiffusion);
    state.(pe).(co).(am).SolidDiffusion = model.(pe).(co).(am).(sd).updateAverageConcentration(state.(pe).(co).(am).SolidDiffusion);
    
    cAverageNE = state.(ne).(co).(am).(sd).cAverage;
    cAveragePE = state.(pe).(co).(am).(sd).cAverage;
    cAverageElyte  = state.(elyte).c;
    
    vol_part_ne = (4/3).* pi .* rp0_ne .^3;
    vol_part_pe = (4/3).* pi .* rp0_pe .^3;
    vols_ne = model.(ne).(co).grid.cells.volumes;
    vols_pe = model.(pe).(co).grid.cells.volumes;
    volsElyte = model.(elyte).grid.cells.volumes;

    Npart_ne = vf_ne / ((4/3) .* pi.* rp0_ne.^3);
    Npart_pe = vf_pe / ((4/3) .* pi.* rp0_pe .^3);


    NLi_ne = sum(Npart_ne .* vols_ne .* cAverageNE .* vol_part_ne);
    NLi_pe = sum(Npart_pe .* vols_pe .* cAveragePE .* vol_part_pe);
    NLi_elyte = sum(volsElyte .* porosElyte .* cAverageElyte);

    state.(ne).(co).(am) = model.(ne).(co).(am).updateRvol(state.(ne).(co).(am));
    state.(pe).(co).(am) = model.(pe).(co).(am).updateRvol(state.(pe).(co).(am));

    Rvol_ne = state.(ne).(co).(am).Rvol;
    Rvol_pe = state.(ne).(co).(am).Rvol;

    Term 
    
    Y(end+1) = NLi_ne;
    X(end+1) = T(t)/hour;
end


plot(X, Y);

ylabel('total lithium quantity')
xlabel('Time (in hours)')


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












