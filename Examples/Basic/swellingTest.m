%gen
%model static (l 550)
%N*n (l 395)
%list grave ?
%
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

fname = fullfile('ParameterData','BatteryCellParameters',...
                 'LithiumIonBatteryCell','lithium_ion_battery_nmc_silicon.json');
jsonstruct = parseBattmoJson(fname);
% jsonstruct2 = parseBattmoJson("sample_input.json");
% jsonstruct = mergeJsonStructs({jsonstruct, jsonstruct2});

%removing temperature effect

jsonstruct.use_thermal = false;

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

%%
jsonstruct.(ne).(co).(am).(sd).N = 5;
jsonstruct.(pe).(co).(am).(sd).N = 5;
%%
%%%
% It is also possible to update the properties of this inputparams in a
% similar way to updating the jsonstruct. Here we set the discretisation
% level for the diffusion model. Other input parameters for the full diffusion
% model can be found here:
% :class:`FullSolidDiffusionModelInputParams <Electrochemistry.FullSolidDiffusionModelInputParams>`

%%% Setting up the geometry
% Here, we setup the 1D computational mesh that will be used for the
% simulation. The required discretization parameters are already included
% in the class BatteryGenerator1D. Classes for generating other geometries can
% be found in the BattMo/Battery/BatteryGeometry folder.

jsonstruct_geometry = parseBattmoJson('Examples/jsondatafiles/geometry1d.json');
jsonstruct = mergeJsonStructs({jsonstruct_geometry , jsonstruct});

% merging with another json to add a TimeStepping field (good way)?
jsonstruct2 = parseBattmoJson("Examples/Documentation/jsonfiles/Example/timeStepping.json");
jsonstruct = mergeJsonStructs({jsonstruct, jsonstruct2});
%we choose a pertinent totalTime (the simulation stops as the battery is 
% empty, 3960 = 1.1*3600 to get total time in seconds)
jsonstruct.TimeStepping.totalTime = 3960/jsonstruct.Control.DRate;
%set a timestep of 100
jsonstruct.TimeStepping.numberOfTimeSteps = jsonstruct.TimeStepping.totalTime/100;
%set up model and output
output = runBatteryJson(jsonstruct, 'includeGridGenerator', true);
model = output.model;
cgt = model.cgt;
gen = output.gridGenerator;

%% Plot Open Circuit Potential as a funcion of State of Charge for pe and ne
%(using the report p12)

Temp = 298.15;
elde = {ne,pe};

figure
hold on
for i = 1:numel(elde)
   po_itf = model.(elde{i}).(co).(am).(itf);

   theta100 = po_itf.guestStoichiometry100;
   theta0   = po_itf.guestStoichiometry0;
   cmax     = po_itf.saturationConcentration;

   soc   = linspace(0, 1);
   theta = soc*theta100 + (1 - soc)*theta0;
   c     = theta.*cmax;
   OCP = po_itf.computeOCPFunc(c, Temp, cmax);

   plot(soc, OCP)
end
xlabel('SOC [-]')
ylabel('OCV [V]')
title('OCV for both electrodes');
legend(elde)

%% Plotting the results

states = output.states;

%removing the empty fields
ind = cellfun(@(x) ~isempty(x), states);
states = states(ind);

%collecting E and I to plot them in the next cells
E = cellfun(@(x) x.Control.E, states); 
I = cellfun(@(x) x.Control.I, states);

T = cellfun(@(x) x.time, states); 
I;
%% Plot E as a function of the time
figure()
subplot(2,2,1)
plot(T/hour, E)
xlabel('time [hours]')
ylabel('Cell Voltage [V]')

%% Plot I as a function of the time (seems to be constant in this case ?)
subplot(2,2,2)
plot(T/hour, I)
xlabel('time [hours]')
ylabel('Cell Current [A]')

%% Plot the overpotential of the negative electrode as a function of time
subplot(2,2,3)
negativeElectrodeSize = model.NegativeElectrode.G.topology.cells.num;
L = "x = 1";
hold on
for i = 1:negativeElectrodeSize
    phi = cellfun(@(state) state.NegativeElectrode.Coating.phi(i), states);
    phiElyte = cellfun(@(x) x.Electrolyte.phi(i), states);

    cSurf = cellfun(@(x) x.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion.cSurface(i), states);
    
    %assuming saturationConcentration is cmax
    cmax = model.NegativeElectrode.Coating.ActiveMaterial.Interface.saturationConcentration;
    Temp = 298.15;
    
    %again !! assuming saturationConcentration is cmax
    OCP   = model.NegativeElectrode.Coating.ActiveMaterial.Interface.computeOCPFunc(cSurf, Temp, cmax);
    eta = phi - phiElyte - OCP;
    plot(T/hour, eta);
    if i > 1
        L(end+1) = "x = " + int2str(i);
    end
end
xlabel('time [hours]')
ylabel('Eta of the Negative electrode')
legend(L);
hold off
%% Plot the porosity as a function of the time for different position across the negative electrode
%constant of x and time. Pretty sus imo
subplot(2,2,4)
negativeElectrodeSize = model.NegativeElectrode.G.topology.cells.num;
L = "x = 1";
hold on;
for i = 1:negativeElectrodeSize
   %found a porosity in model.Separator.porosity, but it doesn't seem to
   %fit for what we want
   porosity = cellfun(@(x) x.NegativeElectrode.Coating.porosity(i), states);
   plot(T/hour, porosity);
   if i > 1
       L(end+1) = "x = " + int2str(i);
   end
end
xlabel('time [hours]')
ylabel('Porosity of the Negative Electrode')
legend(L);


%% Plot the porosity as a function of the position for different times

subplot(2,2,4)
negativeElectrodeSize = model.NegativeElectrode.G.topology.cells.num;
N_elements_ne = gen.nenx;
deltaX = (negativeElectrodeSize/(N_elements_ne-1)) * 10^6;
totalTime = length(T);
initstate = output.states{1}; %seems like an assumption
position = (0:N_elements_ne-1) * deltaX; %tab pos[i] = (i-1)*deltaX

legendTime = "t = 0 hour";
porosity = zeros(1, N_elements_ne);
for i = 1:N_elements_ne
    porosity(i) = initstate.NegativeElectrode.Coating.porosity(i);
end
plot(position, porosity);

hold on
for t = 1:totalTime
    if t < 80
        timestep = 8;
    else
        timestep = 110;
    end

    if mod(t, timestep) == 0
        porosity = initstate.NegativeElectrode.Coating.porosity(1:N_elements_ne);
        plot(position, porosity);
        hold on; 
        if t > 1
            currentTimeHour = T(t) / hour;
            legendTime(end+1) = "t = " + num2str(currentTimeHour, 2) + " hour";
        end
    end
end
hold off
xlabel('Position across the Negative Electrode (in µm)')
ylabel('Porosity of the Negative Electrode')
legend(legendTime);

%% Plot the porosity near the cc as a function of the state of charge
%something is plotted but constant... Gotta see the model
%subplot(2,2,4)
figure
 
 ne_itf = model.NegativeElectrode.(co).(am).(itf);
 ne_sd  = model.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion;
 
 totalTime = length(T);
 realTotalTime = states{totalTime}.time;
 
 theta100 = ne_itf.guestStoichiometry100;
 theta0   = ne_itf.guestStoichiometry0;
 cmax     = ne_itf.saturationConcentration;
 r0       = ne_sd.particleRadius; %actually wait is it ??
 F        = model.con.F;
 N        = ne_sd.N;
 
 
 soc = [];
 porosity = [];
 
 for t = 1:totalTime
     sumConcentrations = 0;
     for i = 1:N
         c = states{t}.NegativeElectrode.(co).ActiveMaterial.SolidDiffusion.c(i);
         sumConcentrations = sumConcentrations + c;
     end
     cAverage = sumConcentrations/N;
 
     theta = cAverage/cmax;
 
     poros = states{t}.NegativeElectrode.Coating.porosity(1); % warning taking the first one. 
     % And I don't know what it means. The vector is constant anyway lmao
     
     stoc = (theta-theta0)./(theta100-theta0);
 
     soc(end+1) = stoc;
     porosity(end+1) = poros;
 end
 porosit
 plot(soc, porosity,'LineWidth',1.2);
 axis([0 1 0 1]);
 axis square;
 
 xlabel('State of Charge near the current collector')
 ylabel('Porosity near the current collector')

%% Plot the concentration of the electrolyte as a function of the position for different times
figure
%subplot(2,2,4)
negativeElectrodeSize = model.NegativeElectrode.G.topology.cells.num;
positiveElectrodeSize = model.PositiveElectrode.G.topology.cells.num;
separatorSize         = model.Separator.G.topology.cells.num;
N_elements_ne = gen.nenx;
N_elements_pe = gen.penx;
N_elements_sep = gen.sepnx;

deltaX_ne  = (negativeElectrodeSize/(N_elements_ne-1)) * 10^6;
deltaX_pe  = (positiveElectrodeSize/(N_elements_pe-1)) * 10^6;
deltaX_sep = (separatorSize/(N_elements_sep-1)) * 10^6;

totalTime = length(T);

%constructing the x axis
%beep boop. Assembling the 3 zones.
pos_ne = (0:N_elements_ne-1) * deltaX_ne;
pos_sep = negativeElectrodeSize * 1e6 + (1:N_elements_sep-1) * deltaX_sep;
pos_pe = (negativeElectrodeSize + separatorSize) * 1e6 + (1:N_elements_pe-1) * deltaX_pe;

position = [pos_ne, pos_sep, pos_pe];


legendTime = "t = 0 hour";

%beep boop again
c_ne  = c(1:N_elements_ne);
c_sep = c(N_elements_ne + (1:N_elements_sep-1));
c_pe  = c(N_elements_ne + N_elements_sep + (1:N_elements_pe-1));

%Assemble
concentration = [c_ne, c_sep, c_pe];


plot(position, concentration);

hold on

for t = 1:totalTime
    disp(totalTime)
   %Only draw the curve for timestep multiples
   if t<75
       timestep = 8;
   else
       timestep = 110;
   end

   if mod(t,timestep) == 0
       %In the initial code, he excluded the upper bound of each part
       N_total = N_elements_ne + N_elements_sep + N_elements_pe;
       concentration = states{t}.Electrolyte.c(1:N_total);
       plot(position, concentration);
       if t > 1
           t = T(t)/hour;
           legendTime(end+1) = "t = " + num2str(t,2) + " hour";
       end

   end 
end
%hold on necessary ??? on the previous code yup
xline(negativeElectrodeSize* 10^6,'red','separator');
xline((negativeElectrodeSize + separatorSize)* 10^6,'red');

xlabel('Position across the Electrolyte (in µm)')
ylabel('Concentration of the Electrolyte (in mol/m3')
legend(legendTime);



%% Plot the radius evolution

figure

ne_itf   = model.NegativeElectrode.(co).(am).(itf);
ne_sd = model.NegativeElectrode.(co).ActiveMaterial.SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

theta100 = ne_itf.guestStoichiometry100;
theta0   = ne_itf.guestStoichiometry0;
cmax     = ne_itf.saturationConcentration;
r0       = ne_sd.particleRadius;
F        = model.con.F;
N        = ne_sd.N; %dim radius - c is 100*1 sphere by sphere
N_elements_ne = 10;  %gen.nenx; %heh


Y = [];
X = [];

for t = 1:totalTime

    sumTheta = 0;
        for x = 1:N_elements_ne       
        sumConcentrations = 0;
        for i = 1:N
            c = states{t}.(ne).(co).(am).(sd).c((x-1)*N +i);
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

ne_itf   = model.NegativeElectrode.Coating.ActiveMaterial.Interface;
ne_sd = model.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

theta100 = ne_itf.guestStoichiometry100;
theta0   = ne_itf.guestStoichiometry0;
cmax     = ne_itf.saturationConcentration;
r0       = ne_sd.particleRadius;
F        = model.con.F;
N        = ne_sd.N;
N_elements_ne = 10; %gen.nenx; %still meh


Y = [];
X = [];

for t = 1:totalTime


    state = states{t};
    state.(ne).(co).(am).(sd) = ne_sd.updateAverageConcentration(state.(ne).(co).(am).(sd));
    cAverage = state.(ne).(co).(am).(sd).cAverage;

    vols = model.NegativeElectrode.Coating.ActiveMaterial.G.cells.volumes;

    theta = (sum(vols .* cAverage))./(sum(vols).* cmax);

    sumTheta = 0;
       for x = 1:N_elements_ne       
       sumConcentrations = 0;
       for i = 1:N
           c = states{t}.NegativeElectrode.(co).ActiveMaterial.SolidDiffusion.c((x-1)*N +i);
           sumConcentrations = sumConcentrations + c;
       end
       cAverage = sumConcentrations/N;

       theta = cAverage/cmax;

       sumTheta = sumTheta + theta;
       end
       theta = sumTheta/N_elements_ne;

       soc = (theta - theta0)/(theta100-theta0);


    Y(end+1) = soc;
    X(end+1) = T(t)/hour;
end

plot(X, Y);

ylabel('State of Charge in the NEGATIVE ELECTRODE')
xlabel('Time (in hours)')

%% now the positive electrode
figure
subplot(2,1,2)

po_itf   = model.PositiveElectrode.(co).(am).(itf);
po_sd = model.PositiveElectrode.(co).ActiveMaterial.SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

theta100 = po_itf.guestStoichiometry100;
theta0   = po_itf.guestStoichiometry0;
cmax     = po_itf.saturationConcentration;
r0       = po_sd.particleRadius; %this variable is not used. Anyway why does it exist??
F        = model.con.F;
N        = po_sd.N;
N_elements_pe = 10; %gen.penx;


Y = [];
X = [];

for t = 1:totalTime

    sumTheta = 0;
    for x = 1:N_elements_pe
        sumConcentrations = 0;
        for i = 1:N
            c = states{t}.PositiveElectrode.Coating.ActiveMaterial.SolidDiffusion.c((x-1)*N + i);
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

ne_itf   = model.NegativeElectrode.Coating.(am).(itf);
ne_sd = model.NegativeElectrode.Coating.ActiveMaterial.SolidDiffusion;

totalTime = length(T);
realTotalTime = states{totalTime}.time;

theta100 = ne_itf.guestStoichiometry100;
theta0   = ne_itf.guestStoichiometry0;
cmax     = ne_itf.saturationConcentration;
r0       = ne_sd.particleRadius;
F        = model.con.F;
N        = ne_sd.N;
N_elements_ne = 10; %gen.nenx; %still meh


Y = [];
X = [];

vf_ne = model.NegativeElectrode.(co).(am).(sd).volumeFraction;
vf_pe = model.PositiveElectrode.(co).(am).(sd).volumeFraction;
porosElyte = model.Electrolyte.volumeFraction;
rp0_ne = model.NegativeElectrode.(co).ActiveMaterial.SolidDiffusion.particleRadius;
rp0_pe = model.PositiveElectrode.(co).ActiveMaterial.SolidDiffusion.particleRadius;

for t = 1:totalTime


    state = states{t};

    state.NegativeElectrode.(co).(am).(sd) = model.(ne).(co).(am).SolidDiffusion.updateAverageConcentration(state.(ne).(co).(am).SolidDiffusion);
    state.PositiveElectrode.(co).(am).(sd) = model.(pe).(co).(am).SolidDiffusion.updateAverageConcentration(state.(pe).(co).(am).SolidDiffusion);
    
    cAverageNE = state.(ne).(co).(am).SolidDiffusion.cAverage;
    cAveragePE = state.(pe).(co).(am).SolidDiffusion.cAverage;
    cAverageElyte  = state.Electrolyte.c;
    
    vol_part_ne = (4/3).* pi .* rp0_ne .^3;
    vol_part_pe = (4/3).* pi .* rp0_pe .^3;
    vols_ne = model.(ne).(co).(am).G.cells.volumes;
    vols_pe = model.(pe).(co).(am).G.cells.volumes;
    volsElyte = model.Electrolyte.G.cells.volumes;

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












