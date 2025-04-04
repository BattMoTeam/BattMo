
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

%%% Specifying the physical model
% In this tutorial we will simulate a lithium-ion battery consisting of a 
% negative electrode, a positive electrode and an electrolyte. *BattMo* 
% comes with some pre-defined models which can be loaded from JSON files.
% Here we will load the basic lithium-ion model JSON file which comes with
% Battmo.

fname = fullfile('ParameterData','BatteryCellParameters',...
                 'LithiumIonBatteryCell','lithium_ion_battery_nmc_silicon.json');
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
% jsonstruct. To make life easier we define some shorthand
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
%%
jsonstruct.(ne).(co).(am).(sd).N = 5;
jsonstruct.(pe).(co).(am).(sd).N = 5;
jsonstruct.(pe).(co).(am).(sd)
%%
%%%
% It is also possible to update the properties of this inputparams in a
% similar way to updating the jsonstruct. Here we set the discretisation
% level for the diffusion model. Other input parameters for the full diffusion
% model can be found here:
% :class:`FullSolidDiffusionModelInputParams <Electrochemistry.FullSolidDiffusionModelInputParams>`.

%%inputparams.(ne).(co).(am).(sd).N = 5;
%%inputparams.(pe).(co).(am).(sd).N = 5;

%%% Setting up the geometry
% Here, we setup the 1D computational mesh that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGenerator1D. Classes for generating other geometries can
% be found in the BattMo/Battery/BatteryGeometry folder.

jsonstruct_geometry = parseBattmoJson('Examples/jsondatafiles/geometry1d.json');
jsonstruct = mergeJsonStructs({jsonstruct_geometry , jsonstruct});

%set up model and output
output = runBatteryJson(jsonstruct);
model = output.model
cgt = model.cgt;
%%Plot Open Circuit Potential as a funcion of State of Charge for pe and ne
%%(using the report p12)
%%%
% Now, we update the inputparams with the properties of the mesh. This function
% will update relevent parameters in the inputparams object and make sure we have
% all the required parameters for the model geometry chosen.
%%inputparams = gen.updateBatteryInputParams(inputparams);
%%% Initialising the battery model object
% The battery model is initialized by sending inputparams to the Battery class
% constructor. see :class:`Battery <Battery.Battery>`.
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

% T = 298.15;
% elde = {ne,pe};
% 
% figure
% hold on
% for i = 1:numel(elde)
%    po_itf = model.(elde{i}).(co).(am).(itf);
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
% end
% xlabel('SOC [-]')
% ylabel('OCV [V]')
% title('OCV for both electrodes');
% legend(elde)

%%% Controlling the simulation
% The control model specifies how the simulation is controlled. This can
% also be thought of as the boundary conditions of the simulation.

%%%
% In the first instance we use IEswitch control policy.
% We set the total time scaled by the CRate in the model.
% The CRate has been set by the json file. We can access it here:
%%
timestep.timeStepDuration = 100;

step    = model.Control.setupScheduleStep(timestep);
control = model.Control.setupScheduleControl();

%This control is used to set up the schedule
%Noote : double switch structure seems obsolete, as initialControl is containent in controlPolicy now
%gotta find the Control Policy possibilities
schedule = struct('control', control, 'step', step);
schedule = model.Control.setupSchedule
%step contient un array de int qui renvoie à control (array)

tup = 0.1;
switch model.Control.controlPolicy
  case {'IEswitch'}
      switch model.Control.initialControl
            case 'discharging'
                inputI = model.Control.Imax;
                inputE = model.Control.lowerCutoffVoltage;
            case {'charging', 'first-charge'}
                inputI = -model.Control.Imax;
                inputE = model.Control.upperCutoffVoltage;
            otherwise
        error('initCase not recognized')
      end
      srcfunc = @(time, I, E) rampupSwitchControl(time, tup, I, E, ...
                                                inputI, ...
                                                inputE);
      control.IEswitch = true;
      control = struct('src', srcfunc, 'IEswitch', true);
  case 'CC'
    srcfunc = @(time) rampupControl(time, tup, model.Control.Imax);
    control.CC = true;
  otherwise
    error('control policity not recognized');
end
%%
%%%
% We create a control structure containing the source function and
% specifying that we want to use IESwitch control:

% control.src = srcfunc;

%%% Setting the initial state of the battery
% To run simulation we need to know the starting point which we will run it
% from, in terms of the value of the primary variables being modelled at
% the start of the simulation. 
% The initial state of the model is setup using model.setupInitialState()
% Here we take the state of charge (SOC) given in the input and calculate
% equilibrium concentration based on theta0, theta100 and cmax.

%Noote : he's trying to set up an initial state, already done with the new
%battmo ? never seen this before
initstate = model.setupInitialState(jsonstruct);


%%% Running the simulation
% Once we have the initial state, the model and the schedule, we can call
% the simulateScheduleAD function which will actually run the simulation:

%model.verbose = true;

nls = NonLinearSolver;
nls.errorOnFailure = false;


[wellSols, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls); 

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
subplot(2,2,3)
negativeElectrodeSize = model.NegativeElectrode.G.cells.num;
L = "x = 1";
for i = 1:negativeElectrodeSize
    hold on

    phi = cellfun(@(x) x.NegativeElectrode.ActiveMaterial.phi(i), states);
    phiElyte = cellfun(@(x) x.Electrolyte.phi(i), states);

    cSurf = cellfun(@(x) x.NegativeElectrode.ActiveMaterial.SolidDiffusion.cSurface(i), states);
    cmax = model.NegativeElectrode.ActiveMaterial.Interface.cmax;
    Temp = 298.15;

    
    OCP   = model.NegativeElectrode.ActiveMaterial.Interface.computeOCPFunc(cSurf, Temp, cmax);
    eta = phi - phiElyte - OCP;
    plot(T/hour, eta);


    if i > 1
        L(end+1) = "x = " + int2str(i);
    end
end
xlabel('time [hours]')
ylabel('Eta of the Negative electrode')
legend(L);

%% Plot the porosity as a function of the time for different position across the negative electrode
%subplot(2,2,4)
%negativeElectrodeSize = model.NegativeElectrode.G.cells.num;
%L = "x = 1";
%for i = 1:negativeElectrodeSize
%    hold on
%    porosity = cellfun(@(x) x.NegativeElectrode.ActiveMaterial.porosity(i), states);
%    plot(T/hour, porosity);
%    if i > 1
%        L(end+1) = "x = " + int2str(i);
%    end
%end
%xlabel('time [hours]')
%ylabel('Porosity of the Negative Electrode')
%legend(L);



% Plot the porosity as a function of the position for different times
subplot(2,2,4)
negativeElectrodeSize = gen.xlength(2);
N_elements_ne = gen.nenx;
deltaX = (negativeElectrodeSize/(N_elements_ne-1)) * 10^6;
totalTime = length(T);

position = [];
for x = 1:N_elements_ne
    position(end+1) = (x-1)*deltaX;
end

legendTime = "t = 0 hour";
porosity = [];
for i = 1:N_elements_ne
            porosity(end+1) = initstate.NegativeElectrode.ActiveMaterial.porosity(i);
end
plot(position, porosity);


for t = 1:totalTime
    hold on
    %Only draw the curve for timestep multiples
    if t<80
        timestep = 8;
    else
        timestep = 110;
    end
    if mod(t,timestep) == 0
        porosity = [];
        for i = 1:N_elements_ne
            porosity(end+1) = states{t}.NegativeElectrode.ActiveMaterial.porosity(i);
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
%
%
%
%
 %% Plot the porosity near the cc as a function of the state of charge
 
%subplot(2,2,4)
figure
 
 po_itf   = model.NegativeElectrode.(am).(itf);
 ne_sd = model.NegativeElectrode.ActiveMaterial.SolidDiffusion;
 
 totalTime = length(T);
 realTotalTime = states{totalTime}.time;
 
 theta100 = po_itf.theta100;
 theta0   = po_itf.theta0;
 cmax     = po_itf.cmax;
 r0       = ne_sd.rp;
 F        = model.con.F;
 N        = ne_sd.N;
 
 
 soc = [];
 porosity = [];
 
 for t = 1:totalTime
     sumConcentrations = 0;
     for i = 1:N
         c = states{t}.NegativeElectrode.ActiveMaterial.SolidDiffusion.c(i);
         sumConcentrations = sumConcentrations + c;
     end
     cAverage = sumConcentrations/N;
 
     theta = cAverage/cmax;
 
     poros = states{t}.NegativeElectrode.ActiveMaterial.porosity(1);
     
     stoc = (theta-theta0)./(theta100-theta0);
 
     soc(end+1) = stoc;
     porosity(end+1) = poros;
 end
 
 plot(soc, porosity,'LineWidth',1.2);
 axis([0 1 0 1]);
 axis square;
 
 xlabel('State of Charge near the current collector')
 ylabel('Porosity near the current collector')
%% Plot the concentration of the electrolyte as a function of the position for different times
%figure
%%subplot(2,2,4)
%negativeElectrodeSize = gen.xlength(2);
%positiveElectrodeSize = gen.xlength(4);
%separatorSize          = gen.xlength(3);
%
%N_elements_ne = gen.nenx;
%N_elements_pe = gen.penx;
%N_elements_sep = gen.sepnx;
%
%deltaX_ne  = (negativeElectrodeSize/(N_elements_ne-1)) * 10^6;
%deltaX_pe  = (positiveElectrodeSize/(N_elements_pe-1)) * 10^6;
%deltaX_sep = (separatorSize/(N_elements_sep-1)) * 10^6;
%
%totalTime = length(T);
%
%%constructing the x axis
%position = [];
%for x = 1:N_elements_ne
%    position(end+1) = (x-1)*deltaX_ne;
%end
%for x = 1:N_elements_sep-1
%    position(end+1) = negativeElectrodeSize* 10^6 + x*deltaX_sep;
%end
%for x = 1:N_elements_pe-1
%    position(end+1) = (negativeElectrodeSize + separatorSize) * 10^6 + x*deltaX_pe;
%end
%
%
%legendTime = "t = 0 hour";
%concentration = [];
%for i = 1:N_elements_ne
%            concentration(end+1) = initstate.Electrolyte.c(i);
%end
%for i = 1:N_elements_sep-1
%            concentration(end+1) = initstate.Electrolyte.c(N_elements_ne + i);
%end
%for i = 1:N_elements_pe-1
%            concentration(end+1) = initstate.Electrolyte.c(N_elements_ne + N_elements_sep + i);
%end
%
%
%plot(position, concentration);
%
%
%for t = 1:totalTime
%    hold on
%    %Only draw the curve for timestep multiples
%    if t<75
%        timestep = 8;
%    else
%        timestep = 110;
%    end
%
%    if mod(t,timestep) == 0
%        concentration = [];
%        for i = 1:N_elements_ne
%            concentration(end+1) = states{t}.Electrolyte.c(i);
%        end
%        for i = 1:N_elements_sep-1
%            concentration(end+1) = states{t}.Electrolyte.c(N_elements_ne + i);
%        end
%        for i = 1:N_elements_pe-1
%            concentration(end+1) = states{t}.Electrolyte.c(N_elements_ne + N_elements_sep + i);
%        end
%        
%        plot(position, concentration);
%
%        if t > 1
%            t = T(t)/hour;
%            legendTime(end+1) = "t = " + num2str(t,2) + " hour";
%        end
%
%    end 
%end
%
%hold on
%xline(negativeElectrodeSize* 10^6,'red','separator');
%xline((negativeElectrodeSize + separatorSize)* 10^6,'red');
%
%xlabel('Position across the Electrolyte (in µm)')
%ylabel('Concentration of the Electrolyte (in mol/m3')
%legend(legendTime);



%% Plot the radius evolution

figure

ne_itf   = model.NegativeElectrode.(am).(itf);
ne_sd = model.NegativeElectrode.ActiveMaterial.SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

theta100 = ne_itf.theta100;
theta0   = ne_itf.theta0;
cmax     = ne_itf.cmax;
r0       = ne_sd.rp;
F        = model.con.F;
N        = ne_sd.N;
N_elements_ne = gen.nenx;


Y = [];
X = [];

for t = 1:totalTime

    sumTheta = 0;
        for x = 1:N_elements_ne       
        sumConcentrations = 0;
        for i = 1:N
            c = states{t}.NegativeElectrode.ActiveMaterial.SolidDiffusion.c((x-1)*N +i);
            sumConcentrations = sumConcentrations + c;
        end
        cAverage = sumConcentrations/N;
    
        radius = computeRadius(cAverage,cmax,r0);
        end
    
    Y(end+1) = radius;
    X(end+1) = T(t)/hour;
end

plot(X, Y);
ylabel('Silicon particle radius')
xlabel('Time (in hours)')












%% Plot soc for each electrode as a function of time

figure
subplot(2,1,1)

ne_itf   = model.NegativeElectrode.(am).(itf);
ne_sd = model.NegativeElectrode.ActiveMaterial.SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

theta100 = ne_itf.theta100;
theta0   = ne_itf.theta0;
cmax     = ne_itf.cmax;
r0       = ne_sd.rp;
F        = model.con.F;
N        = ne_sd.N;
N_elements_ne = gen.nenx;


Y = [];
X = [];

for t = 1:totalTime


    state = states{t};
    state.NegativeElectrode.ActiveMaterial.SolidDiffusion = model.NegativeElectrode.ActiveMaterial.SolidDiffusion.updateAverageConcentration(state.NegativeElectrode.ActiveMaterial.SolidDiffusion);
    cAverage = state.NegativeElectrode.ActiveMaterial.SolidDiffusion.cAverage;

    vols = model.NegativeElectrode.ActiveMaterial.G.cells.volumes;

    theta = (sum(vols .* cAverage))./(sum(vols).* cmax);

    %sumTheta = 0;
    %    for x = 1:N_elements_ne       
    %    sumConcentrations = 0;
    %    for i = 1:N
    %        c = states{t}.NegativeElectrode.ActiveMaterial.SolidDiffusion.c((x-1)*N +i);
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
        soc = (theta - theta0)/(theta100-theta0);


    Y(end+1) = soc;
    X(end+1) = T(t)/hour;
end

plot(X, Y);

ylabel('State of Charge in the NEGATIVE ELECTRODE')
xlabel('Time (in hours)')



subplot(2,1,2)

po_itf   = model.PositiveElectrode.(am).(itf);
po_sd = model.PositiveElectrode.ActiveMaterial.SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

theta100 = po_itf.theta100;
theta0   = po_itf.theta0;
cmax     = po_itf.cmax;
r0       = po_sd.rp;
F        = model.con.F;
N        = po_sd.N;
N_elements_pe = gen.penx;


Y = [];
X = [];

for t = 1:totalTime

    sumTheta = 0;
    for x = 1:N_elements_pe
        sumConcentrations = 0;
        for i = 1:N
            c = states{t}.PositiveElectrode.ActiveMaterial.SolidDiffusion.c((x-1)*N + i);
            sumConcentrations = sumConcentrations + c;
        end
        cAverage = sumConcentrations/N;    
        theta = cAverage/cmax;

        sumTheta = sumTheta + theta;
    end
    theta = sumTheta/N_elements_pe;

    soc = (theta-theta100)/(theta0-theta100);

    Y(end+1) = soc;
    X(end+1) = T(t)/hour;
end

plot(X, Y);

ylabel('State of charge in the POSITIVE ELECTRODE')
xlabel('Time (in hours)')






%%

figure

ne_itf   = model.NegativeElectrode.(am).(itf);
ne_sd = model.NegativeElectrode.ActiveMaterial.SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

theta100 = ne_itf.theta100;
theta0   = ne_itf.theta0;
cmax     = ne_itf.cmax;
r0       = ne_sd.rp;
F        = model.con.F;
N        = ne_sd.N;
N_elements_ne = gen.nenx;


Y = [];
X = [];

for t = 1:totalTime

    vf_ne = model.NegativeElectrode.ActiveMaterial.Interface.volumeFraction;
    vf_pe = model.PositiveElectrode.ActiveMaterial.Interface.volumeFraction;
    porosElyte = model.Electrolyte.volumeFraction;
    rp0_ne = model.NegativeElectrode.ActiveMaterial.SolidDiffusion.rp;
    rp0_pe = model.PositiveElectrode.ActiveMaterial.SolidDiffusion.rp;



    state = states{t};

    state.NegativeElectrode.ActiveMaterial.SolidDiffusion = model.NegativeElectrode.ActiveMaterial.SolidDiffusion.updateAverageConcentration(state.NegativeElectrode.ActiveMaterial.SolidDiffusion);
    state.PositiveElectrode.ActiveMaterial.SolidDiffusion = model.PositiveElectrode.ActiveMaterial.SolidDiffusion.updateAverageConcentration(state.PositiveElectrode.ActiveMaterial.SolidDiffusion);
    
    cAverageNE = state.NegativeElectrode.ActiveMaterial.SolidDiffusion.cAverage;
    cAveragePE = state.PositiveElectrode.ActiveMaterial.SolidDiffusion.cAverage;
    cAverageElyte  = state.Electrolyte.c;
    
    vol_part_ne = (4/3).* pi .* rp0_ne .^3;
    vol_part_pe = (4/3).* pi .* rp0_pe .^3;
    vols_ne = model.NegativeElectrode.ActiveMaterial.G.cells.volumes;
    vols_pe = model.PositiveElectrode.ActiveMaterial.G.cells.volumes;
    volsElyte = model.Electrolyte.G.cells.volumes;

    Npart_ne = vf_ne / ((4/3) .* pi.* rp0_ne.^3);
    Npart_pe = vf_pe / ((4/3) .* pi.* rp0_pe .^3);


    NLi_ne = sum(Npart_ne .* vols_ne .* cAverageNE .* vol_part_ne);
    NLi_pe = sum(Npart_pe .* vols_pe .* cAveragePE .* vol_part_pe);
    NLi_elyte = sum(volsElyte .* porosElyte .* cAverageElyte);



    state.NegativeElectrode.ActiveMaterial = model.NegativeElectrode.ActiveMaterial.updateRvol(state.NegativeElectrode.ActiveMaterial);
    state.PositiveElectrode.ActiveMaterial = model.PositiveElectrode.ActiveMaterial.updateRvol(state.PositiveElectrode.ActiveMaterial);

    Rvol_ne = state.NegativeElectrode.ActiveMaterial.Rvol;
    Rvol_pe = state.NegativeElectrode.ActiveMaterial.Rvol;

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












