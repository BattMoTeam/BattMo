% 1D test

clear all
close all

%% add MRST module
mrstModule add ad-core

% setup object and run simulation
% delete('temp2.txt');
obj = lithiumIonModel1D();

%% run simulation

[t, y] = obj.p2d();

%% plotting

mrstModule add mrst-gui

fv = obj.fv;

%%  plot of potential

for iy = 1 : size(y, 1)
    yy = y(iy, :)';
    varname = 'E';
    E(iy) = yy(fv.getSlot(varname));
end

%% 
figure
plot(t, E)
title('Potential (E)')
xlabel('time')


