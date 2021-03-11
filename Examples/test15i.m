% _2D_test_case

clear all
close all

%% Add MRST module
mrstModule add ad-core multimodel mrst-gui battery
mrstVerbose off

modelcase = '2D';
switch modelcase
  case '2D'
    inputparams = BatteryInputParams2D();
  case '3D'
    inputparams = BatteryInputParams3D();
end

model = BatteryModel15i(inputparams);

%% run simulation

[t, y] = model.runSimulation(); 

%% process output data

initstate = model.setupInitialState(); 
model = setupFV(model, initstate); 
sl = model.fv.slots; 

clear state E
for iy = 1 : size(y, 1)
    yy = y(iy, :)'; 
    states_elyte{iy}.Li = yy(sl{1}); 
    states_elyte{iy}.phi = yy(sl{2}); 
    states_ne{iy}.Li = yy(sl{3}); 
    states_ne{iy}.phi = yy(sl{4}); 
    states_pe{iy}.Li = yy(sl{5}); 
    states_pe{iy}.phi = yy(sl{6}); 
    states_ccne{iy}.phi = yy(sl{7}); 
    states_ccpe{iy}.phi = yy(sl{8}); 
    E(iy) = yy(sl{9}); 
end

plot(t/hour, E)
