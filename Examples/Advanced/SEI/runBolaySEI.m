%% Particle simulation with SEI layer growth (Bolay et al 2022)

% clear the workspace and close open figures
clear
close all

%% Import the required modules from MRST
% load MRST modules
mrstModule add ad-core mrst-gui mpfa

ne    = 'NegativeElectrode';
pe    = 'PositiveElectrode';
am    = 'ActiveMaterial';
co    = 'Coating';
sd    = 'SolidDiffusion';
itf   = 'Interface';
sei   = 'SolidElectrodeInterface';
sr    = 'SideReaction';
elyte = 'Electrolyte';
ctrl  = 'Control';

%% Setup the properties of the Li-ion battery materials and of the cell design
jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                        'lithium_ion_battery_nmc_graphite.json');
jsonstruct = parseBattmoJson(jsonfilename);

jsonstruct.use_thermal = false;

jsonfilename = fullfile('ParameterData', 'ParameterSets', 'Bolay2022', 'bolay_sei_interface.json');
jsonstruct_bolay = parseBattmoJson(jsonfilename);


jsonstruct.(ne).(co).(am) = mergeJsonStructs({jsonstruct.(ne).(co).(am), ...
                                              jsonstruct_bolay});

jsonstruct.(ne).(co).(am).SEImodel = 'Bolay';

jsontruct_control = struct( 'controlPolicy'     , 'CCCV'       , ...
                            'initialControl'    , 'discharging', ...
                            'numberOfCycles'    , 1            , ...
                            'CRate'             , 1            , ...
                            'DRate'             , 1            , ...
                            'lowerCutoffVoltage', 3            , ...
                            'upperCutoffVoltage', 4            , ...
                            'dIdtLimit'         , 1e-2         , ...
                            'dEdtLimit'         , 1e-2);

jsonstruct.(ctrl) = jsontruct_control;

%% 

doCompare = false;

if doCompare
    
    jsonfilename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', ...
                            'lithium_ion_battery_nmc_graphite.json');
    jsonstruct_original = parseBattmoJson(jsonfilename);

    jsonstruct_original.use_thermal = false;

    jsonfilename = fullfile('ParameterData', 'ParameterSets', 'Bolay2022', 'bolay_sei_interface.json');
    jsonstruct_original_bolay = parseBattmoJson(jsonfilename);


    jsonstruct_original.(ne).(co).(am) = mergeJsonStructs({jsonstruct_original.(ne).(co).(am), ...
                                                  jsonstruct_original_bolay});

    jsonstruct_original.(ne).(co).(am).SEImodel = 'Bolay';

    jsonstruct_original.(ctrl) = jsontruct_control;


    fjv = compareJson(jsonstruct, jsonstruct_original);

    fjv = fjv.sort('comparison')
    
    fjv.print();
    
    return
    
end
%%

doPrintJsonStruct = false;

if doPrintJsonStruct
    fjv = flattenJsonStruct(jsonstruct);
    % fjv.print('filter', {'parame name', 'SEI'})
end

%%

inputparams = BatteryInputParams(jsonstruct);

gen = BatteryGeneratorP2D();
inputparams = gen.updateBatteryInputParams(inputparams);

model = GenericBattery(inputparams);

%% Setup the schedule
%

% schedule = model.(ctrl).setupSchedule([]);
jsonstruct.TimeStepping.timeStepDuration = 200;
schedule = model.(ctrl).setupSchedule(jsonstruct);

%% Setup the initial state of the model
% The initial state of the model is setup using the model.setupInitialState() method.

initstate = model.setupInitialState();

%% Setup non-linear solver

nls = NonLinearSolver();
nls.errorOnFailure = false;

model.nonlinearTolerance = 1e-5;

%% Run simulation

model.verbose = true;
[~, states, report] = simulateScheduleAD(initstate, model, schedule, 'OutputMinisteps', true, 'NonLinearSolver', nls);

%% Setup for plotting

ind = cellfun(@(x) not(isempty(x)), states); 
states = states(ind);
time = cellfun(@(x) x.time, states); 
E    = cellfun(@(x) x.Control.E, states); 
I    = cellfun(@(x) x.Control.I, states);

for istate = 1 : numel(states)
    states{istate} = model.addVariables(states{istate});
end

%% Plotting

close all

set(0, 'defaultlinelinewidth', 3)
set(0, 'defaultaxesfontsize', 15)

%%

figure
plot(time/hour, E, '*-');
title('Voltage / V')
xlabel('Time / h')

%%

figure
plot(time/hour, I);
title('Current / A')
xlabel('Time / h')

%%

figure
hold on

delta = cellfun(@(state) state.(ne).(co).(am).(itf).SEIlength(end), states);
plot(time/hour, delta/(nano*meter), 'displayname', 'at x_{max}')

delta = cellfun(@(state) state.(ne).(co).(am).(itf).SEIlength(1), states);
plot(time/hour, delta/(nano*meter), 'displayname', 'at x_{min}')

title('SEI thickness in negative electrode/ nm')
xlabel('Time / h')

legend show

%%

figure
hold on

u = cellfun(@(state) state.(ne).(co).(am).(itf).SEIvoltageDrop(end), states);
plot(time/hour, u, 'displayname', 'at x_{max}')

u = cellfun(@(state) state.(ne).(co).(am).(itf).SEIvoltageDrop(1), states);
plot(time/hour, u, 'displayname', 'at x_{min}')

title('SEI voltage drop in negative electrode/ V')
xlabel('Time / h')

%%
figure

vols = model.(ne).(co).G.getVolumes();

for istate = 1 : numel(states)
    state = states{istate};
    m(istate) = sum(vols.*state.(ne).(co).(am).(sd).cAverage);
end

plot(time/hour, m);
title('total lithium amount in negative electrode / mol')
xlabel('Time / h')

legend show

%%

quantities = [];

for timeindex = 1 : numel(states)

	l = states{timeindex}.(ne).(co).(am).(itf).SEIlength;

        indAm = model.(ne).(co).compInds.(am);
        
	vols = model.(ne).(co).G.getVolumes();

        vsa = model.(ne).(co).(am).(itf).volumetricSurfaceArea;
        
	scoef   = model.(ne).(co).(am).(itf).SEIstochiometricCoeffcient;
	seimvol = model.(ne).(co).(am).(itf).SEImolarVolume;

        Liqqt = (scoef/seimvol)*sum(l.*(vsa*vols));
	quantities(end + 1) = Liqqt;
        
end

figure 
plot(time/hour, quantities);
title('Lithium quantity consummed');
xlabel('Time / h');
ylabel('quantity / mol');
grid on;

%%

PE_Li_quantities          = [];
NE_Li_quantities          = [];
Electrolyte_Li_quantities = [];
Electrodes_Li_quantities  = [];
Total_Li_quantities       = [];

for timeindex = 1 : numel(states)
    
	PE_qtt = sum(model.(pe).(co).volumeFractions(1).*model.(pe).(co).volumeFraction.*model.(pe).G.getVolumes.*states{timeindex}.(pe).(co).(am).(sd).cAverage);
	NE_qtt = sum(model.(ne).(co).volumeFractions(1).*model.(ne).(co).volumeFraction.*model.(ne).G.getVolumes.*states{timeindex}.(ne).(co).(am).(sd).cAverage);
	Elyte_qtt = sum(model.Electrolyte.volumeFraction.*model.Electrolyte.G.getVolumes.*states{timeindex}.Electrolyte.c);

	Elode_qtt = PE_qtt + NE_qtt;
	Tot_Liqqt = PE_qtt + NE_qtt + Elyte_qtt;
	
	PE_Li_quantities(end + 1)          = PE_qtt;
	NE_Li_quantities(end + 1)          = NE_qtt;
	Electrolyte_Li_quantities(end + 1) = Elyte_qtt;
	Electrodes_Li_quantities(end + 1)  = Elode_qtt;
	Total_Li_quantities(end + 1)       = Tot_Liqqt;

end

%%

figure
hold on

plot(time/hour, PE_Li_quantities,'DisplayName','Positive Electrode');
plot(time/hour, NE_Li_quantities,'DisplayName','Negative Electrode');
plot(time/hour, Electrolyte_Li_quantities,'DisplayName','Electrolyte');
plot(time/hour, Electrodes_Li_quantities,'DisplayName','Both Electrodes');
plot(time/hour, Total_Li_quantities,'DisplayName','Total (except SEI)');
plot(time/hour, Total_Li_quantities + quantities,'DisplayName','Total (including SEI)');
plot(time/hour, quantities,'DisplayName','In the SEI');
title('Lithium quantity');
xlabel('Time / h');
ylabel('quantity / mol');
grid on;
legend show

%%

figure 
plot(time/hour, 100 * quantities / Total_Li_quantities(1));
title('Percentage of Lithium consummed');
xlabel('Time / h');
ylabel('%');
grid on;
