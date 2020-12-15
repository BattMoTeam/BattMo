clear all
close all

%% add MRST module
mrstModule add ad-core

% setup object and run simulation
% delete('temp2.txt');
obj = lithiumIonModel1D();

%% run simulation

[t, y] = obj.p2d();

%%

% obj = lithiumIon();
% [t, y] = obj.p2d();



