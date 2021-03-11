% _2D_test_case

clear all
close all

%% Add MRST module
mrstModule add ad-core multimodel mrst-gui
mrstVerbose off

model = BatteryModel3Dextra15i();
% model = BatteryModel();
% model.verbose = true;
model.J = 0.1; 



%% run simulation
% profile - detail builtin

[t, y] = model.p2d(); 

initstate = icp2d(model); 
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
