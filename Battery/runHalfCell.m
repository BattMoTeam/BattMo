jsonstruct = parseBattmoJson(fullfile('ParameterData','BatteryCellParameters','LithiumIonBatteryCell','lithium_ion_battery_nmc_graphite.json'));

% We define some shorthand names for simplicity.
elyte   = 'Electrolyte';
am      = 'ActiveMaterial';
itf     = 'Interface';
sd      = 'SolidDiffusion';
ctrl    = 'Control';

jsonstruct.use_thermal = false;
jsonstruct.include_current_collectors = false;

jsonstruct.(pe).(am).diffusionModelType = 'full';
jsonstruct.(ne).(am).diffusionModelType = 'full';

paramobj = BatteryInputParams(jsonstruct);

gen = BatteryGenerator1D();

% Now, we update the paramobj with the properties of the mesh. 
paramobj = gen.updateBatteryInputParams(paramobj);
G = paramobj.(ne).(am).G;

batterymodel = Battery(paramobj);
batterymodel = batterymodel.setupComputationalGraph();
cgtbattery = batterymodel.computationalGraph;

hcjsonstruct.(am) = jsonstruct.(ne).(am);
hcjsonstruct.(elyte).cInit   = 1*mol/litre;
hcjsonstruct.(elyte).phiInit = 0;

hcjsonstruct.(ctrl) = jsonstruct.(ctrl);


paramobj = HalfCellInputParams(hcjsonstruct);
paramobj.(am).G = G;

model = HalfCell(paramobj);

model = model.setupComputationalGraph();
cgt = model.computationalGraph;
[g, edgelabels] = getComputationalGraph(cgt);
plot(g)


% setup initstate

% setup schedule


[~, states, report] = simulateScheduleAD(initstate, model, schedule); 


